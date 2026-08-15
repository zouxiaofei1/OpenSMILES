"""Benchmark preview API: run predictions against a selectable data file.

Data source: data/<data_file>.json (default merged_benchmark.json, 4070 gold rows).
Cache:       tools/benchmark_pred_preview_<stem>.json (legacy name kept for default).

GET  /api/v1/benchmark-preview/datasets   — list benchmark-shaped data files under data/
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

from server import history_store

router = APIRouter(prefix="/api/v1", tags=["benchmark"])

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
DEFAULT_DATA_FILE = "merged_benchmark.json"
CACHE = ROOT / "tools" / "benchmark_pred_preview_data.json"

# Subprocess handle for background generation. Only one live generation runs at
# a time; _proc_data records which data file it targets.
_proc: subprocess.Popen | None = None
_proc_data: str | None = None
_proc_lock = threading.Lock()
_gen_total = 0
_captured: list[str] = []  # last stdout lines for diagnostics


def _resolve_source(data_file: str | None) -> Path:
    """Resolve a benchmark data file under data/, with traversal guard.

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


def _source_total(source: Path) -> int:
    """Return the number of rows in the data file, or 0."""
    if not source.is_file():
        return 0
    try:
        rows = json.loads(source.read_text(encoding="utf-8"))
        return len(rows) if isinstance(rows, list) else 0
    except Exception:
        return 0


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
    src_total = _source_total(_resolve_source(data_file))
    job = history_store.get_job(full, _hist_kind(data_file))
    generating = job is not None and job.proc is not None and job.proc.poll() is None
    return {
        "rows": rows,
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
    out: list[str] = []
    for p in sorted(DATA_DIR.glob("*.json")):
        if not p.is_file():
            continue
        try:
            obj = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(obj, list) or not obj:
            continue
        first = obj[0]
        if isinstance(first, dict) and "smiles" in first and "english_name" in first:
            out.append(p.name)
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

    return {
        "rows": rows,
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

    script = f'''
import json, sys, time, os
from pathlib import Path
sys.path.insert(0, r"{wt / 'src'}")

# Suppress RDKit C++ warnings — they flood stderr and block the pipe buffer
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh
SOURCE = Path(r"{source}")
CACHE = Path(r"{cache_path}")
_SIG = f"{{SOURCE.stat().st_size}}:{{SOURCE.stat().st_mtime_ns}}"

def score_pred(pe, pz, ge, gz):
    en_ok = bool(ge) and normalize_en(pe) == normalize_en(ge)
    zh_ok = bool(gz) and normalize_zh(pz) == normalize_zh(gz)
    if ge and gz: dual = en_ok and zh_ok
    elif ge: dual = en_ok
    elif gz: dual = zh_ok
    else: dual = False
    return {{"en_ok": en_ok, "zh_ok": zh_ok, "ok": dual, "ret": bool(pe or pz)}}

rows = json.loads(SOURCE.read_text(encoding="utf-8"))
total = len(rows)
namer = SMILESNNamer()
payload = []
t0 = time.perf_counter()
for i, row in enumerate(rows):
    smi = str(row.get("smiles") or "")
    ge = row.get("english_name") or ""
    gz = row.get("chinese_name") or ""
    try:
        r = namer.name(smi)
        pe, pz = r.en or "", r.zh or ""
    except Exception:
        pe, pz = "", ""
    sc = score_pred(pe, pz, ge, gz)
    payload.append({{"s": smi, "en": pe, "zh": pz, "ge": ge, "gz": gz, **sc}})
    if (i + 1) % 200 == 0 or (i + 1) == total:
        tmp = CACHE.with_suffix(".tmp")
        tmp.write_text(json.dumps({{"_data_sig": _SIG, "rows": payload}}, ensure_ascii=False), encoding="utf-8")
        tmp.replace(CACHE)
        elapsed = time.perf_counter() - t0
        print(f"PROGRESS {{i+1}}/{{total}} {{(i+1)/elapsed:.1f}}r/s", flush=True)

elapsed = time.perf_counter() - t0
n_ok = sum(1 for r in payload if r["ok"])
print(f"DONE {{total}} rows {{elapsed:.1f}}s dual_ok={{n_ok}}", flush=True)
'''

    try:
        proc = subprocess.Popen(
            [sys.executable, "-c", script],
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

    # Build a small inline script that does the generation
    script = f'''
import json, sys, time
from pathlib import Path
sys.path.insert(0, r"{ROOT / 'src'}")

# Suppress RDKit C++ warnings — they flood stderr and block the pipe buffer
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh
SOURCE = Path(r"{source}")
CACHE = Path(r"{cache_path}")

def score_pred(pe, pz, ge, gz):
    en_ok = bool(ge) and normalize_en(pe) == normalize_en(ge)
    zh_ok = bool(gz) and normalize_zh(pz) == normalize_zh(gz)
    if ge and gz: dual = en_ok and zh_ok
    elif ge: dual = en_ok
    elif gz: dual = zh_ok
    else: dual = False
    return {{"en_ok": en_ok, "zh_ok": zh_ok, "ok": dual, "ret": bool(pe or pz)}}

rows = json.loads(SOURCE.read_text(encoding="utf-8"))
total = len(rows)
namer = SMILESNNamer()
payload = []
t0 = time.perf_counter()
for i, row in enumerate(rows):
    smi = str(row.get("smiles") or "")
    ge = row.get("english_name") or ""
    gz = row.get("chinese_name") or ""
    try:
        r = namer.name(smi)
        pe, pz = r.en or "", r.zh or ""
    except Exception:
        pe, pz = "", ""
    sc = score_pred(pe, pz, ge, gz)
    payload.append({{"s": smi, "en": pe, "zh": pz, "ge": ge, "gz": gz, **sc}})
    if (i + 1) % 200 == 0 or (i + 1) == total:
        tmp = CACHE.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
        tmp.replace(CACHE)
        elapsed = time.perf_counter() - t0
        print(f"PROGRESS {{i+1}}/{{total}} {{(i+1)/elapsed:.1f}}r/s", flush=True)

elapsed = time.perf_counter() - t0
n_ok = sum(1 for r in payload if r["ok"])
print(f"DONE {{total}} rows {{elapsed:.1f}}s dual_ok={{n_ok}}", flush=True)
'''

    try:
        _proc = subprocess.Popen(
            [sys.executable, "-c", script],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,  # RDKit warnings go to stderr; suppress to avoid pipe blocking
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
