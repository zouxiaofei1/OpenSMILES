"""Call graph API: dynamic cProfile sampling of the namer pipeline.

GET /api/v1/call-graph — runs cProfile over the first N benchmark molecules in a
subprocess (so disk code changes are always reflected), parses the pstats profile,
keeps only the requested module's functions, and returns the call graph as
nodes/edges with timing weights for the frontend force-directed graph.

Why a subprocess: the project is installed editable (pyproject `where=["src"]`), so
the server process may hold already-imported modules that would not refresh when
disk source changes. A fresh `python -c` subprocess re-imports current code, making
the mtime-based cache invalidation meaningful.

The subprocess script embeds logic ported from tools/profile_callchain.py (sampling)
and tools/render_callchain.py (pstats filtering + w normalization).
"""

from __future__ import annotations

import colorsys
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException, Query

from server import history_store

router = APIRouter(prefix="/api/v1", tags=["call-graph"])

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "merged_benchmark.json"

_CACHE: dict[str, dict[str, Any]] = {}
_CACHE_LOCK = threading.Lock()
_MAX_CACHE_KEYS = 5
_SRC_SIG_CACHE: tuple[list[tuple[str, int, int]], str] | None = None

# graphviz 分层 SVG 渲染（dot.exe）
_DOT_EXE_CANDIDATES = (r"C:\Program Files\Graphviz\bin\dot.exe",)
_MAX_SVG_NODES = 2500
_SVG_CACHE: dict[str, str] = {}
_MAX_SVG_KEYS = 10
_DOT_EXE: str | None = None

# live sampling progress for the call-graph UI progress bar
_PROGRESS_PREFIX = "CALLGRAPH_PROGRESS "
_PROGRESS: dict[str, Any] = {"done": 0, "total": 0, "active": False}
_PROGRESS_LOCK = threading.Lock()


def _set_progress(done: int, total: int, active: bool) -> None:
    with _PROGRESS_LOCK:
        _PROGRESS.update(done=done, total=total, active=active)

def _source_total() -> int:
    """Number of rows in the benchmark data file, or 0."""
    if not DATA.is_file():
        return 0
    try:
        rows = json.loads(DATA.read_text(encoding="utf-8"))
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

def _src_signature() -> str:
    """sha256 over (relpath, size, mtime_ns) of all src/namepredict/**/*.py.

    File list is cached; when the list is unchanged the hash is reused.
    """
    global _SRC_SIG_CACHE
    base = ROOT / "src" / "namepredict"
    files: list[tuple[str, int, int]] = []
    for f in sorted(base.rglob("*.py")):
        try:
            st = f.stat()
        except OSError:
            continue
        files.append((f.relative_to(ROOT).as_posix(), st.st_size, st.st_mtime_ns))
    if _SRC_SIG_CACHE is not None and _SRC_SIG_CACHE[0] == files:
        return _SRC_SIG_CACHE[1]
    h = hashlib.sha256()
    for item in files:
        h.update(repr(item).encode("utf-8"))
    sig = h.hexdigest()[:16]
    _SRC_SIG_CACHE = (files, sig)
    return sig

def _data_signature() -> str:
    try:
        st = DATA.stat()
        return f"{st.st_size}:{st.st_mtime_ns}"
    except OSError:
        return "missing"

def _get_raw(
    n: int, module: str, aggregate_external: bool, refresh: bool, commit: str | None = None,
) -> tuple[dict[str, Any], bool]:
    """Fetch raw sampled graph from cache, sampling in a subprocess on miss.

    Returns (raw, cached). raw may be {"ok": False, "error": ...} on failure.
    Shared by /call-graph (JSON) and /call-graph/svg so both see the same data.
    With a commit, results are cached under tools/history_cache and keyed by
    sampling params (commit code is immutable, so no src signature needed).
    """
    if commit is not None:
        return _get_raw_history(n, module, aggregate_external, refresh, commit)

    key = (
        f"n={n};mod={module};agg={aggregate_external};"
        f"src={_src_signature()};data={_data_signature()}"
    )
    with _CACHE_LOCK:
        raw = None
        if not refresh:
            raw = _CACHE.get(key)
        if raw is None:
            _set_progress(0, 0, True)
            raw = _run_profile_subprocess(n, module, aggregate_external)
            if not raw.get("ok"):
                _set_progress(0, 0, False)
                return raw, False
            raw["generated_at"] = datetime.now(timezone.utc).isoformat()
            _CACHE[key] = raw
            while len(_CACHE) > _MAX_CACHE_KEYS:
                _CACHE.pop(next(iter(_CACHE)))
            return raw, False
        return raw, True


def _get_raw_history(
    n: int, module: str, aggregate_external: bool, refresh: bool, commit: str,
) -> tuple[dict[str, Any], bool]:
    """Fetch raw graph for a past commit from its history cache (sample on miss)."""
    key = f"n={n};mod={module};agg={aggregate_external}"
    cache = history_store.read_cache(commit, "callgraph")
    if cache and cache.get("_data_sig") == history_store.data_sig():
        entries = cache.get("entries") or {}
        if not refresh and key in entries:
            raw = dict(entries[key])
            raw["ok"] = True
            return raw, True
    _set_progress(0, 0, True)
    raw = _run_profile_subprocess(n, module, aggregate_external, commit)
    if not raw.get("ok"):
        _set_progress(0, 0, False)
        return raw, False
    raw["generated_at"] = datetime.now(timezone.utc).isoformat()
    cache = history_store.read_cache(commit, "callgraph") or {
        "_data_sig": history_store.data_sig(), "entries": {},
    }
    entries = cache.setdefault("entries", {})
    entries[key] = raw
    history_store.write_cache(commit, "callgraph", cache)
    return raw, False

_PROFILE_SAMPLER = ROOT / "server" / "profile_sampler.py"

def _run_profile_subprocess(
    n: int, module: str, aggregate_external: bool, commit: str | None = None,
) -> dict[str, Any]:
    """Run the multi-process sampler script; parse its JSON graph from stdout.

    Reads stdout line-by-line so the sampler's CALLGRAPH_PROGRESS lines update
    the global _PROGRESS dict that /call-graph/progress reports to the UI.

    With a commit: the sampler runs against that commit's worktree src (via
    NAMEPREDICT_SRC_ROOT), the worktree is ref-counted (released in finish_job),
    and the timeout is relaxed to 300s (historical sampling is slower).
    """
    cmd = [sys.executable, str(_PROFILE_SAMPLER), "--n", str(n), "--module", module]
    if aggregate_external:
        cmd.append("--agg")
    env: dict[str, str] | None = None
    timeout = 30
    if commit is not None:
        try:
            wt = history_store.ensure_worktree(commit)  # refs+1; released in finish_job
        except RuntimeError as exc:
            return {"ok": False, "error": f"worktree: {exc}"}
        if not (wt / "src" / "namepredict").is_dir():
            # The package is an editable install pointing at the main repo, so an
            # import would silently fall back to current code — reject instead.
            history_store.release_worktree(commit)
            return {"ok": False, "error": f"该 commit 无 namepredict 源码（{commit[:7]}），无法采样历史调用链"}
        cmd += ["--data", str(DATA)]
        env = dict(os.environ)
        env["NAMEPREDICT_SRC_ROOT"] = str(wt / "src")
        timeout = 300
    lines: list[str] = []
    try:
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            cwd=str(ROOT),
            env=env,
        )
    except OSError as exc:
        if commit is not None:
            history_store.release_worktree(commit)
        return {"ok": False, "error": f"profile subprocess failed to start: {exc}"}
    if commit is not None:
        ok, err = history_store.register_job(commit, "callgraph", p, n)
        if not ok:
            try:
                p.kill()
            except Exception:
                pass
            history_store.release_worktree(commit)
            return {"ok": False, "error": err}
    for line in p.stdout:
        if line.startswith(_PROGRESS_PREFIX):
            parts = line.strip().split()
            if len(parts) == 3:
                try:
                    _set_progress(int(parts[1]), int(parts[2]), True)
                except ValueError:
                    pass
            continue
        lines.append(line)
    try:
        p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        if commit is not None:
            history_store.finish_job(commit, "callgraph")
        return {"ok": False, "error": "profile timed out"}
    if commit is not None:
        history_store.finish_job(commit, "callgraph")
    _set_progress(_PROGRESS["done"], _PROGRESS["total"], False)
    out = "".join(lines)
    marker = "__CALLGRAPH_JSON__"
    if marker not in out:
        tail = out[-500:]
        return {"ok": False, "error": f"profile subprocess failed: {tail}"}
    try:
        raw = json.loads(out.split(marker, 1)[1])
    except Exception as exc:
        return {"ok": False, "error": f"bad profile output: {exc}"}
    raw["ok"] = True
    return raw

def _filter_and_build(
    raw: dict[str, Any], floor_pct: float
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Re-index ids after dropping nodes below floor_pct (edges crossing dropped
    nodes are removed). w normalization is kept global so colors don't jump with
    the threshold."""
    keep = [nd for nd in raw["nodes"] if nd["cum_pct"] >= floor_pct]
    idmap = {old: i for i, old in enumerate(nd["id"] for nd in keep)}
    nodes = []
    for nd in keep:
        n = dict(nd)
        n["id"] = idmap[nd["id"]]
        nodes.append(n)
    edges = []
    for e in raw["edges"]:
        if e["from"] in idmap and e["to"] in idmap:
            e2 = dict(e)
            e2["id"] = len(edges)
            e2["from"] = idmap[e["from"]]
            e2["to"] = idmap[e["to"]]
            edges.append(e2)
    return nodes, edges

def _find_dot() -> str | None:
    global _DOT_EXE
    if _DOT_EXE:
        return _DOT_EXE
    for cand in _DOT_EXE_CANDIDATES:
        if Path(cand).is_file():
            _DOT_EXE = cand
            return cand
    _DOT_EXE = shutil.which("dot")
    return _DOT_EXE

def _heat(w: float) -> str:
    """Blue→red heat color (ported from tools/render_callchain.py heat())."""
    h = (1.0 - min(1.0, max(0.0, w))) * 0.66
    r, g, b = colorsys.hsv_to_rgb(h, 0.85, 0.92)
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def _hash_bits(s: str) -> int:
    """确定性字符串哈希（md5 前 4 字节），跨进程/机器颜色稳定。"""
    return int.from_bytes(hashlib.md5(s.encode()).digest()[:4], "big")


def _dir_hue(module: str) -> float:
    """节点色相：按完整 module 路径哈希，在 0-360° 均匀分布。

    早期按 layerN 主色分段、目录作层内偏移；但实际 module 多形如
    layer2/arene_carbonyl.py（文件直接挂 layer 下、无子目录），目录解析后只剩
    layer2，同层文件色相完全一致、图上糊成一片。改为全 module 哈希后，每个
    文件/文件夹在色相环上都有独立且均匀散布的颜色，层内 20+ 文件也可分辨。
    """
    m = module or "?"
    return (_hash_bits(m) % 360) / 360.0


def _node_color(module: str, w: float) -> str:
    """文件色系节点色：色相由 module 决定，耗时 w 越大越亮/越饱和。"""
    hue = _dir_hue(module)
    light = 0.42 + 0.38 * max(0.0, min(1.0, w))
    sat = 0.62 + 0.30 * max(0.0, min(1.0, w))
    r, g, b = colorsys.hls_to_rgb(hue, light, sat)
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def _text_color(hex_color: str) -> str:
    """按背景亮度选白/黑字，保证可读。"""
    try:
        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)
    except (ValueError, IndexError):
        return "#ffffff"
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return "#ffffff" if lum < 150 else "#0f172a"

def _dot_escape(s: str) -> str:
    """Escape backslashes/quotes for dot string literals (paths/names)."""
    return s.replace("\\", "\\\\").replace('"', '\\"')

def _short_node_name(nd: dict) -> str:
    mod = (nd.get("module") or "").rsplit("/", 1)[-1].replace(".py", "")
    line = nd.get("line")
    label = nd.get("label") or ""
    return f"{mod}:{line}:{label}" if line is not None else f"{mod}:{label}"

def _nodes_edges_to_dot(nodes: list[dict], edges: list[dict]) -> str:
    """Build a graphviz dot from filtered nodes/edges (ids already re-indexed).

    Nodes with external=True (aggregated cross-layer calls) are drawn grey/dashed
    with a (L{n}) tag; cross-layer edges (external=True) are dashed.
    """
    lines = [
        "digraph {",
        '\tgraph [fontname=Arial, nodesep=0.125, ranksep=0.25, bgcolor="transparent"];',
        '\tnode [fontcolor=white, fontname=Arial, height=0, shape=box, style=filled, width=0];',
        '\tedge [fontname=Arial, fontcolor="#334155"];',
    ]
    for nd in nodes:
        tip = _dot_escape(nd.get("file") or nd.get("module") or "")
        if nd.get("external"):
            tag = f"(L{nd.get('layer')})" if nd.get("layer") is not None else "(外部)"
            label = (
                f"{_dot_escape(tag)} {_dot_escape(nd['label'] or '')}\\n"
                f"{nd['cum_pct']:.1f}%\\n{nd['ncalls']}×"
            )
            lines.append(
                f'\t{nd["id"]} [id="cg{nd["id"]}", color="#94a3b8", fontcolor="#1e293b", '
                f'style="filled,dashed", fontsize="10.00", label="{label}", tooltip="{tip}"];'
            )
        else:
            name = _dot_escape(_short_node_name(nd))
            label = (
                f"{name}\\n"
                f"{nd['cum_pct']:.1f}%\\n"
                f"({nd['self_pct']:.2f}%)\\n{nd['ncalls']}×"
            )
            bg = _node_color(nd.get("module") or "", nd.get("w") or 0.0)
            fg = _text_color(bg)
            lines.append(
                f'\t{nd["id"]} [id="cg{nd["id"]}", color="{bg}", fontcolor="{fg}", '
                f'fontsize="11.00", label="{label}", tooltip="{tip}"];'
            )
    for e in edges:
        w = max(0.5, e["w"] * 5)
        label = f'{e["cum_pct"]:.1f}%\\n{e["calls"]}×'
        style = ', style="dashed"' if e.get("external") else ""
        lines.append(
            f'\t{e["from"]} -> {e["to"]} [label="{label}", penwidth="{w:.2f}"{style}];'
        )
    lines.append("}")
    return "\n".join(lines)

def _layer_subgraph(
    nodes: list[dict], edges: list[dict], sel
) -> tuple[list[dict], list[dict]]:
    """裁剪出选中 layer 的子图，跨层调用聚合为 external 节点（虚线灰）。

    sel: callable(node) -> bool 判断节点是否属于目标 layer。输出节点的 id
    重新编号；external 边带 external=True 标记，外部节点带 external=True。
    """
    inside = [n for n in nodes if sel(n)]
    if not inside:
        return [], []
    inside_ids = {n["id"] for n in inside}
    ext: dict[int, dict] = {}
    merged: list[dict] = []
    for e in edges:
        fi, ti = e["from"], e["to"]
        if fi in inside_ids and ti in inside_ids:
            merged.append(e)
        elif fi in inside_ids or ti in inside_ids:
            e2 = dict(e)
            e2["external"] = True
            merged.append(e2)
            for nid in (fi, ti):
                if nid not in inside_ids and nid not in ext:
                    orig = next((n for n in nodes if n["id"] == nid), None)
                    if orig:
                        ext[nid] = orig
    new_id: dict[int, int] = {}
    new_nodes: list[dict] = []
    for n in inside:
        new_id[n["id"]] = len(new_nodes)
        nn = dict(n)
        nn["id"] = len(new_nodes)
        nn["external"] = False
        new_nodes.append(nn)
    for oid, on in ext.items():
        new_id[oid] = len(new_nodes)
        nn = dict(on)
        nn["id"] = len(new_nodes)
        nn["external"] = True
        new_nodes.append(nn)
    new_edges: list[dict] = []
    for e in merged:
        if e["from"] in new_id and e["to"] in new_id:
            ne = dict(e)
            ne["id"] = len(new_edges)
            ne["from"] = new_id[e["from"]]
            ne["to"] = new_id[e["to"]]
            new_edges.append(ne)
    return new_nodes, new_edges

def _render_svg(dot_text: str) -> tuple[str | None, str | None]:
    dot_exe = _find_dot()
    if not dot_exe:
        return None, "graphviz dot 未找到（需安装或配置 PATH）"
    try:
        p = subprocess.run(
            [dot_exe, "-Tsvg"],
            input=dot_text,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:
        return None, f"dot 渲染失败: {exc}"
    if p.returncode != 0:
        return None, f"dot 退出码 {p.returncode}: {(p.stderr or '')[-500:]}"
    svg = p.stdout
    start = svg.find("<svg")
    end = svg.rfind("</svg>")
    if start == -1 or end == -1:
        return None, "dot 输出缺少 <svg>"
    return svg[start:end + 6], None

@router.get("/call-graph/svg")
def call_graph_svg(
    n: int = Query(200, ge=1, le=100000),
    floor_pct: float = Query(1.0, ge=0, le=100),
    refresh: bool = False,
    module: str = "namepredict",
    aggregate_external: bool = False,
    layer: int | None = Query(None, ge=-1, le=6),
    commit: str | None = Query(None),
) -> dict[str, Any]:
    """Layered graphviz SVG of the sampled call graph (pan/zoom on the frontend).

    layer: None=全部；0..5=仅该 pipeline layer；6=Tools 共享层；
    -1=核心调度（无 layer 归属的函数，如 namer.py 入口）。
    """
    full = _hist_commit(commit)
    if not DATA.is_file():
        return {"ok": False, "error": "benchmark data not found"}
    src_total = _source_total()
    if src_total <= 0:
        return {"ok": False, "error": "benchmark data empty or unreadable"}
    n = min(n, src_total)

    t0 = time.perf_counter()
    raw, cached = _get_raw(n, module, aggregate_external, refresh, full)
    if not raw.get("ok"):
        return raw
    nodes, edges = _filter_and_build(raw, floor_pct)
    if layer is not None:
        if layer == -1:
            nodes, edges = _layer_subgraph(nodes, edges, lambda nd: nd.get("layer") is None)
        else:
            nodes, edges = _layer_subgraph(nodes, edges, lambda nd: nd.get("layer") == layer)
    if len(nodes) > _MAX_SVG_NODES:
        return {
            "ok": False,
            "error": f"节点数 {len(nodes)} 超过上限 {_MAX_SVG_NODES}，"
            f"请提高累计耗时阈值（当前 floor_pct={floor_pct}%）",
        }

    key = (
        f"n={n};mod={module};agg={aggregate_external};floor={floor_pct:.3f};"
        f"layer={layer};src={full[:10] if full else _src_signature()};data={_data_signature()}"
    )
    svg = None
    svg_cached = False
    with _CACHE_LOCK:
        if not refresh:
            svg = _SVG_CACHE.get(key)
            svg_cached = svg is not None
        if svg is None:
            dot_text = _nodes_edges_to_dot(nodes, edges)
            svg, err = _render_svg(dot_text)
            if err:
                return {"ok": False, "error": err}
            _SVG_CACHE[key] = svg
            while len(_SVG_CACHE) > _MAX_SVG_KEYS:
                _SVG_CACHE.pop(next(iter(_SVG_CACHE)))

    meta = {
        "n": n,
        "source_total": src_total,
        "n_actual": raw.get("n_actual", n),
        "total_s": round(raw.get("wall_s", raw.get("total_s", 0.0)), 6),
        "n_nodes_total": len(raw["nodes"]),
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "floor_pct": floor_pct,
        "cached": cached,
        "svg_cached": svg_cached,
        "elapsed_ms": round((time.perf_counter() - t0) * 1000.0, 1),
        "generated_at": raw.get("generated_at", ""),
    }
    return {"ok": True, "svg": svg, "meta": meta}


@router.get("/call-graph/progress")
def call_graph_progress() -> dict[str, Any]:
    """Live sampling progress for the call-graph UI progress bar."""
    with _PROGRESS_LOCK:
        done, total, active = (
            _PROGRESS["done"], _PROGRESS["total"], _PROGRESS["active"],
        )
    pct = round(100.0 * done / total, 1) if total > 0 else 0.0
    return {"ok": True, "done": done, "total": total, "pct": pct, "active": active}


@router.get("/call-graph")
def call_graph(
    n: int = Query(200, ge=1, le=100000),
    floor_pct: float = Query(0.1, ge=0, le=100),
    refresh: bool = False,
    module: str = "namepredict",
    aggregate_external: bool = False,
    commit: str | None = Query(None),
) -> dict[str, Any]:
    """Dynamically sampled call graph of the naming pipeline."""
    full = _hist_commit(commit)
    if not DATA.is_file():
        return {"ok": False, "error": "benchmark data not found"}
    src_total = _source_total()
    if src_total <= 0:
        return {"ok": False, "error": "benchmark data empty or unreadable"}
    n = min(n, src_total)

    t0 = time.perf_counter()
    raw, cached = _get_raw(n, module, aggregate_external, refresh, full)
    if not raw.get("ok"):
        return raw
    nodes, edges = _filter_and_build(raw, floor_pct)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    meta = {
        "n": n,
        "source_total": src_total,
        "n_actual": raw.get("n_actual", n),
        "total_s": round(raw.get("wall_s", raw.get("total_s", 0.0)), 6),
        "n_nodes_total": len(raw["nodes"]),
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "floor_pct": floor_pct,
        "cached": cached,
        "elapsed_ms": round(elapsed_ms, 1),
        "generated_at": raw.get("generated_at", ""),
    }
    return {"ok": True, "nodes": nodes, "edges": edges, "meta": meta}
