"""Benchmark preview API: run predictions against a selectable data file.

Data source: benchmarks/<data_file>.json (default merged_benchmark.json, 4070 gold rows).
Cache:       tools/benchmark_pred_preview_<stem>.json (legacy name kept for default).

GET  /api/v1/benchmark-preview/datasets   — list benchmark-shaped data files under benchmarks/
GET  /api/v1/benchmark-preview        — serve cached rows + generation status
POST /api/v1/benchmark-preview/refresh — start subprocess generation
GET  /api/v1/benchmark-preview/status  — poll generation progress
"""

from __future__ import annotations

import json
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from benchmarks import preview_metrics as metrics
from server.backend import history_store

router = APIRouter(prefix="/api/v1", tags=["benchmark"])

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "benchmarks"
DEFAULT_DATA_FILE = "merged_benchmark.json"
CACHE = ROOT / "tools" / "benchmark_pred_preview_data.json"

# Subprocess handle for background generation. Only one live generation runs at
# a time; _proc_data records which data file it targets.
_proc: subprocess.Popen | None = None
_proc_data: str | None = None
_proc_lock = threading.Lock()
_gen_total = 0
_captured: list[str] = []  # last stdout lines for diagnostics

# 数据集文件索引缓存: (stat 签名, 文件名数组, SMILES → 文件下标数组)。
# 命中的代价只是一次 stat, 所以可以放在每次 GET /benchmark-preview 上。
_idx_lock = threading.Lock()
_idx_state: tuple[tuple, list[str], dict[str, list[int]]] | None = None

# 源文件行缓存: path → (size, mtime_ns, rows)。/status 每 2s 就调一次 _source_total,
# 没有这层缓存时每次都要全量 json.loads 一遍源文件。
_src_lock = threading.Lock()
_src_cache: dict[str, tuple[int, int, list[dict]]] = {}


def _resolve_source(data_file: str | None) -> Path:
    """Resolve a benchmark data file under benchmarks/, with traversal guard.

    Accepts a bare filename (resolved against DATA_DIR). Slashes are rejected
    so callers cannot reach arbitrary files elsewhere on disk.
    """
    name = (data_file or DEFAULT_DATA_FILE).strip()
    if "/" in name or "\\" in name or name in ("", ".", ".."):
        raise ValueError("invalid data file name")
    p = (DATA_DIR / name).resolve()
    if p.suffix.lower() != ".json" or not str(p).startswith(str(DATA_DIR.resolve())):
        raise ValueError("invalid data file path")
    return p


def _cache_path_for(data_file: str) -> Path:
    """Per-file preview cache under tools/, keeping the legacy name for default."""
    if data_file == DEFAULT_DATA_FILE:
        return CACHE
    return ROOT / "tools" / f"benchmark_pred_preview_{Path(data_file).stem}.json"


def _hist_kind(data_file: str) -> str:
    """history_store job/cache kind; default file keeps the legacy key."""
    if data_file == DEFAULT_DATA_FILE:
        return "benchmark"
    return f"benchmark_{Path(data_file).stem}"


def _data_sig_for(source: Path) -> str:
    """Signature of a benchmark data file; cache invalidates when it changes."""
    try:
        st = source.stat()
        return f"{st.st_size}:{st.st_mtime_ns}"
    except OSError:
        return "missing"


def _read_cache(cache: Path) -> list[dict]:
    """Read cached rows, or empty list."""
    if not cache.is_file():
        return []
    try:
        rows = json.loads(cache.read_text(encoding="utf-8"))
        return rows if isinstance(rows, list) else []
    except Exception:
        return []


def _ensure_derived(rows: list[dict], cache: Path) -> list[dict]:
    """为旧缓存 rows 补算派生字段(稠环/复杂度/相似度)并写回, 避免整体重新生成。

    仅当 rows 完整(非生成中)且缺字段时调用; ~5s/4k 行(多数花在 BertzCT 上),
    写回后后续请求直接命中。"""
    if metrics.has_derived(rows):
        return rows
    out_rows = []
    for r in rows:
        out = dict(r)
        out.update(metrics.row_derived_of(out))
        # 排除用的 facet 字段只在响应里存在, 不能跟着派生字段落盘。
        out.pop("feat", None)
        out.pop("f", None)
        out_rows.append(out)
    try:
        tmp = cache.with_suffix(".tmp")
        tmp.write_text(json.dumps(out_rows, ensure_ascii=False), encoding="utf-8")
        tmp.replace(cache)
    except OSError:
        pass  # 写回失败不影响本次返回
    return out_rows


def _stat_key() -> tuple:
    """benchmarks/*.json 的 (name, size, mtime_ns) 快照, 用作索引/源缓存失效判据。

    只 stat 不解析: 新增/删除/改动任一文件都会改变本 key, 因此可以放在请求路径上。
    """
    out = []
    for p in sorted(DATA_DIR.glob("*.json")):
        try:
            st = p.stat()
        except OSError:
            continue
        out.append((p.name, st.st_size, st.st_mtime_ns))
    return tuple(out)


def _is_benchmark_file(p: Path) -> bool:
    """该 json 是否为 benchmark 形状: 非空 list, 且首行是含 smiles/english_name 的 dict。"""
    try:
        obj = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return False
    if not isinstance(obj, list) or not obj:
        return False
    first = obj[0]
    return isinstance(first, dict) and "smiles" in first and "english_name" in first


def _dataset_files() -> list[Path]:
    """benchmarks/ 下所有 benchmark 形状的数据文件, 按文件名排序。"""
    return [p for p in sorted(DATA_DIR.glob("*.json")) if p.is_file() and _is_benchmark_file(p)]


def _file_index() -> tuple[list[str], dict[str, list[int]]]:
    """(数据集文件名数组, SMILES → 所属文件下标数组), 按 stat 签名缓存。

    前端"按文件排除"靠它把行归属回数据集。多文件归属是常态(merged 的行同时属于若干子集),
    所以值恒为下标数组。返回的两半必须来自同一次调用 —— 混用别处算出的文件名数组会让
    下标指向错误的文件。
    """
    global _idx_state
    key = _stat_key()
    with _idx_lock:
        if _idx_state is not None and _idx_state[0] == key:
            return _idx_state[1], _idx_state[2]
        files: list[str] = []
        idx: dict[str, list[int]] = {}
        for p in _dataset_files():
            try:
                obj = json.loads(p.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(obj, list):
                continue
            i = len(files)
            files.append(p.name)
            for r in obj:
                if not isinstance(r, dict):
                    continue
                smi = r.get("smiles")
                if not isinstance(smi, str) or not smi:
                    continue
                bucket = idx.get(smi)
                if bucket is None:
                    idx[smi] = [i]
                elif bucket[-1] != i:
                    # 同一文件内的重复 SMILES 会连续出现, 所以比末项即可去重。
                    bucket.append(i)
        _idx_state = (key, files, idx)
        return files, idx


def _source_rows(source: Path) -> list[dict]:
    """源 benchmark 文件的行列表; 按 (size, mtime_ns) 缓存, 缺失/损坏时返回空表。"""
    try:
        st = source.stat()
    except OSError:
        return []
    key = str(source)
    with _src_lock:
        hit = _src_cache.get(key)
        if hit is not None and hit[0] == st.st_size and hit[1] == st.st_mtime_ns:
            return hit[2]
    try:
        rows = json.loads(source.read_text(encoding="utf-8"))
    except Exception:
        return []
    if not isinstance(rows, list):
        return []
    with _src_lock:
        _src_cache[key] = (st.st_size, st.st_mtime_ns, rows)
    return rows


def _source_total(source: Path) -> int:
    """Return the number of rows in the data file, or 0."""
    return len(_source_rows(source))


def _enrich_rows(rows: list[dict], source: Path, idx: dict[str, list[int]]) -> list[dict]:
    """给预览行补 feat(源行 features) 与 f(所属数据集文件下标), 供前端左栏排除。

    行与源文件行按下标一一对应(见 benchmark_preview_parallel 的 ordering contract), 但
    live 缓存没有 _data_sig 保护: 源文件被插入/重排/等长替换时会整体错位。这时宁可不给
    feat(该维度少过滤), 也不能把别的分子的特征挂到这一行上(会给出错误结果)。首尾行足以
    识别插入/重排/截断这几类错位。

    f 由行自身 SMILES 查表, 与对齐与否无关, 恒可用。rows 来自本次请求刚解析的缓存, 就地
    改字段是安全的(history_store.read_cache / _read_cache 都不共享对象)。
    """
    src = _source_rows(source)
    n = len(src)
    aligned = n >= len(rows) and (
        not rows
        or (rows[0].get("s") == src[0].get("smiles")
            and rows[-1].get("s") == src[len(rows) - 1].get("smiles"))
    )
    for i, r in enumerate(rows):
        feat = src[i].get("features") if aligned and i < n else None
        r["feat"] = [str(x) for x in feat] if isinstance(feat, list) else []
        r["f"] = idx.get(str(r.get("s") or ""), [])
    return rows


def _hist_commit(commit: str | None) -> str | None:
    """Resolve a commit ref → full hash; None = absent/HEAD. 400 on invalid."""
    if not commit:
        return None
    try:
        return history_store.resolve_commit(commit)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"invalid commit: {commit}")


def _hist_rows(full: str, data_file: str) -> list[dict]:
    cache = history_store.read_cache(full, _hist_kind(data_file))
    if cache and cache.get("_data_sig") == _data_sig_for(_resolve_source(data_file)):
        rows = cache.get("rows")
        return rows if isinstance(rows, list) else []
    return []


def _hist_benchmark_preview(full: str, data_file: str) -> dict[str, Any]:
    rows = _hist_rows(full, data_file)
    # 历史缓存可能缺派生字段: 内存补算(不改写历史缓存), 供稠环过滤与相似度排序。
    if not metrics.has_derived(rows):
        rows = [dict(r, **metrics.row_derived_of(r)) for r in rows]
    source = _resolve_source(data_file)
    src_total = _source_total(source)
    files, idx = _file_index()
    rows = _enrich_rows(rows, source, idx)
    job = history_store.get_job(full, _hist_kind(data_file))
    generating = job is not None and job.proc is not None and job.proc.poll() is None
    return {
        "rows": rows,
        "files": files,
        "cached": len(rows),
        "source_total": src_total,
        "stale": len(rows) < src_total,
        "generating": generating,
        "gen_done": job.done if job else len(rows),
        "gen_total": job.total if job else src_total,
        "commit": full,
        "head": history_store._head(),
        "data_file": data_file,
    }


@router.get("/benchmark-preview/datasets")
def benchmark_preview_datasets() -> dict[str, Any]:
    """List data files that look like benchmark rows (list of {smiles, english_name})."""
    out = [p.name for p in _dataset_files()]
    return {"ok": True, "datasets": out}


@router.get("/benchmark-preview")
def get_benchmark_preview(commit: str | None = None, data_file: str | None = None) -> dict[str, Any]:
    """Return benchmark prediction rows from cache + generation status."""
    try:
        source = _resolve_source(data_file)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    data_name = source.name
    full = _hist_commit(commit)
    if full is not None:
        return _hist_benchmark_preview(full, data_name)

    cache = _cache_path_for(data_name)
    rows = _read_cache(cache)
    src_total = _source_total(source)
    generating = _proc is not None and _proc.poll() is None and _proc_data == data_name
    gen_done = len(rows)
    gen_total = _gen_total if _proc_data == data_name else src_total

    # If the subprocess just finished, do a final read
    if not generating and _proc is not None and _proc.poll() == 0 and _proc_data == data_name:
        rows = _read_cache(cache)
        gen_done = len(rows)

    # 旧缓存补算派生字段(一次性写回, 供稠环过滤/相似度排序/复杂度排序)
    if not generating and rows:
        rows = _ensure_derived(rows, cache)

    # facet 必须在 _ensure_derived 写回之后补, 否则会跟着落盘。
    files, idx = _file_index()
    rows = _enrich_rows(rows, source, idx)

    return {
        "rows": rows,
        "files": files,
        "cached": len(rows),
        "source_total": src_total,
        "stale": len(rows) < src_total,
        "generating": generating,
        "gen_done": gen_done,
        "gen_total": gen_total,
        "data_file": data_name,
    }


def _hist_benchmark_refresh(full: str, force: bool, data_file: str) -> dict[str, Any]:
    """Start background generation of benchmark rows for a past commit."""
    source = _resolve_source(data_file)
    src_total = _source_total(source)
    if src_total == 0:
        return {"ok": False, "error": f"source data not found or empty: {source}"}
    try:
        wt = history_store.ensure_worktree(full)  # refs+1; released in finish_job
    except RuntimeError as exc:
        return {"ok": False, "error": str(exc)}
    if not (wt / "src" / "namepredict").is_dir():
        # The package is an editable install pointing at the main repo, so an
        # import would silently fall back to current code — reject instead.
        history_store.release_worktree(full)
        return {"ok": False, "error": f"该 commit 无 namepredict 源码（{full[:7]}），无法运行历史 benchmark"}

    kind = _hist_kind(data_file)
    cache_path = history_store.cache_path(full, kind)
    history_store.cache_dir(full).mkdir(parents=True, exist_ok=True)

    # Parallel generator subprocess (benchmarks.benchmark_preview_parallel):
    # same module as the live preview, but imports namepredict from the commit's
    # worktree (--src) and wraps cache rows with the data signature (--sig).
    cmd = [
        sys.executable, "-m", "benchmarks.benchmark_preview_parallel",
        "--data", str(source),
        "--cache", str(cache_path),
        "--src", str(wt / "src"),
        "--sig",
    ]
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as exc:
        history_store.release_worktree(full)
        return {"ok": False, "error": str(exc)}

    ok, err = history_store.register_job(full, kind, proc, src_total)
    if not ok:
        try:
            proc.kill()
        except Exception:
            pass
        history_store.release_worktree(full)
        return {"ok": False, "error": err, "total": src_total}

    job = history_store.get_job(full, kind)
    _streams_done = [0]

    def _mark_stream_done() -> None:
        _streams_done[0] += 1
        if _streams_done[0] >= 2:
            error = None
            if proc.poll() not in (0, None):
                tail = "".join(job.captured[-10:]) if job and job.captured else ""
                error = tail[-500:] if tail else f"exit code {proc.poll()}"
            history_store.finish_job(full, kind, error=error)

    def _drain_stdout() -> None:
        try:
            if proc.stdout is not None:
                for line in proc.stdout:
                    if job is not None:
                        job.captured.append(line)
                    if line.startswith("PROGRESS"):
                        parts = line.split()
                        if len(parts) >= 2 and "/" in parts[1]:
                            try:
                                if job is not None:
                                    job.done = int(parts[1].split("/")[0])
                            except ValueError:
                                pass
        except Exception:
            pass
        finally:
            _mark_stream_done()

    def _drain_stderr() -> None:
        try:
            if proc.stderr is not None:
                for line in proc.stderr:
                    if job is not None:
                        job.captured.append(line)
        except Exception:
            pass
        finally:
            _mark_stream_done()

    # Separate stdout/stderr drain threads: the historical code can flood stderr
    # with RDKit C++ warnings, which would otherwise fill the pipe buffer and
    # stall the child (the live path avoids this with stderr=DEVNULL).
    threading.Thread(target=_drain_stdout, daemon=True).start()
    threading.Thread(target=_drain_stderr, daemon=True).start()
    return {"ok": True, "started": True, "total": src_total, "commit": full, "data_file": data_file}


@router.post("/benchmark-preview/refresh")
def refresh_benchmark_preview(force: bool = False, commit: str | None = None, data_file: str | None = None) -> dict[str, Any]:
    """Start background subprocess to generate benchmark cache.

    Set force=true to kill any stuck generation and restart.
    """
    try:
        source = _resolve_source(data_file)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    data_name = source.name
    full = _hist_commit(commit)
    if full is not None:
        return _hist_benchmark_refresh(full, force, data_name)

    global _proc, _proc_data, _gen_total, _captured

    # Check for a running process
    if _proc is not None:
        if _proc.poll() is None:
            if force:
                # Kill the stuck process and clean up
                try:
                    _proc.kill()
                except Exception:
                    pass
                _proc = None
                _proc_data = None
                _captured = []
            else:
                hint = "" if _proc_data == data_name else f"（正在生成 {_proc_data}）"
                return {
                    "ok": False,
                    "error": f"generation already in progress{hint} (use force=true to reset)",
                    "total": _gen_total,
                }
        else:
            # Process already dead — clean up the stale reference
            _proc = None
            _proc_data = None
            _captured = []

    src_total = _source_total(source)
    if src_total == 0:
        return {"ok": False, "error": f"source data not found or empty: {source}"}

    _gen_total = src_total
    _captured = []
    _proc_data = data_name
    cache_path = _cache_path_for(data_name)

    # Parallel generator subprocess (benchmarks.benchmark_preview_parallel):
    # process pool over gold rows; keeps cache a continuous prefix + prints
    # PROGRESS n/total so /status and the drain below track progress.
    cmd = [
        sys.executable, "-m", "benchmarks.benchmark_preview_parallel",
        "--data", str(source),
        "--cache", str(cache_path),
    ]
    try:
        _proc = subprocess.Popen(
            cmd,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,  # worker chatter is silenced in the module
            text=True,
        )
    except Exception as exc:
        _proc = None
        _proc_data = None
        return {"ok": False, "error": str(exc)}

    # Drain stdout in a daemon thread so the pipe buffer never blocks the child
    def _drain() -> None:
        if _proc is None or _proc.stdout is None:
            return
        try:
            for line in _proc.stdout:
                _captured.append(line)
        except Exception:
            pass

    threading.Thread(target=_drain, daemon=True).start()

    return {"ok": True, "started": True, "total": src_total, "data_file": data_name}


def _hist_benchmark_status(full: str, data_file: str) -> dict[str, Any]:
    """Poll historical generation progress (same shape as live status)."""
    kind = _hist_kind(data_file)
    done = len(_hist_rows(full, data_file))
    total = _source_total(_resolve_source(data_file))
    job = history_store.get_job(full, kind)

    proc_running = job is not None and job.proc is not None and job.proc.poll() is None
    proc_exited = job is not None and job.proc is not None and job.proc.poll() is not None

    running = proc_running and done < total

    error = None
    if proc_exited:
        if job and job.proc and job.proc.poll() != 0:
            tail = "".join(job.captured[-10:]) if job.captured else ""
            error = tail[-500:] if tail else f"exit code {job.proc.poll()}"
        elif done < total:
            error = "process exited successfully but cache is incomplete"
    if error is None:
        # job removed after finish_job → fall back to the persisted error marker
        error = history_store.read_error(full, kind)

    stuck = not proc_running and not proc_exited and done < total and done > 0 and error is None
    return {"running": running, "done": done, "total": total, "error": error, "stuck": stuck}


@router.get("/benchmark-preview/status")
def benchmark_status(commit: str | None = None, data_file: str | None = None) -> dict[str, Any]:
    """Poll generation progress."""
    try:
        source = _resolve_source(data_file)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    data_name = source.name
    full = _hist_commit(commit)
    if full is not None:
        return _hist_benchmark_status(full, data_name)

    global _proc
    cache = _cache_path_for(data_name)
    done = len(_read_cache(cache))
    total = _gen_total if _proc_data == data_name else _source_total(source)

    proc_running = _proc is not None and _proc.poll() is None and _proc_data == data_name
    proc_exited = _proc is not None and _proc.poll() is not None and _proc_data == data_name

    # If the process exited but cache is incomplete, it failed
    running = proc_running and done < total

    error = None
    if proc_exited:
        if _proc.poll() != 0:
            tail = "".join(_captured[-10:]) if _captured else ""
            error = tail[-500:] if tail else f"exit code {_proc.poll()}"
        # If process exited successfully but done < total, something went wrong
        elif done < total:
            error = "process exited successfully but cache is incomplete"

    # If no process running and cache incomplete, generation is stale/stuck
    stuck = not proc_running and not proc_exited and done < total and done > 0

    return {
        "running": running,
        "done": done,
        "total": total,
        "error": error,
        "stuck": stuck,
    }
