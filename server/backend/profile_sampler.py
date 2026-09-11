"""Multi-process cProfile sampler for the call-graph API.

Why multi-process: cProfile's profiler is thread-local and rdkit is CPU-bound
(GIL), so threading would not parallelise. Each worker profiles a chunk of
molecules, dumps its pstats to a temp file, and the parent merges them.

Windows spawn safety: worker functions (_init_worker/_profile_chunk) are
module-level and the pool is created inside main(), guarded by __main__.

Output: prints "__CALLGRAPH_JSON__" marker then a JSON object
{nodes, edges, total_s, n_actual} — same contract as the old inline script.
"""
from __future__ import annotations

import argparse
import cProfile
import json
import os
import re
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pstats

ROOT = Path(__file__).resolve().parents[1]
# NAMEPREDICT_SRC_ROOT lets the history feature sample a past commit's code from
# its worktree without running that commit's own (possibly missing) sampler.
_SRC_ROOT = Path(os.environ.get("NAMEPREDICT_SRC_ROOT") or (ROOT / "src"))
sys.path.insert(0, str(_SRC_ROOT))

_WORKER_NAMER = None


def _init_worker() -> None:
    global _WORKER_NAMER
    from rdkit import RDLogger

    RDLogger.logger().setLevel(RDLogger.ERROR)
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _report(done: int, total: int) -> None:
    """Emit sampling progress to stdout for the parent process (Popen line reader)."""
    if total <= 0:
        return
    print(f"CALLGRAPH_PROGRESS {done} {total}", flush=True)


def _profile_chunk(smiles: list[str], progress=None) -> str:
    """Profile a chunk of molecules in the current worker; return temp .prof path.

    progress(done, total) is called after each molecule. Only the parent process
    passes a callback (single-process path); workers pass None so their stdout
    stays quiet and the multi-process path reports by completed chunk instead.
    """
    pr = cProfile.Profile()
    pr.enable()
    total = len(smiles)
    for i, s in enumerate(smiles, 1):
        try:
            _WORKER_NAMER.name(s)
        except Exception:
            pass
        if progress is not None:
            progress(i, total)
    pr.disable()
    fd, path = tempfile.mkstemp(suffix=".prof")
    os.close(fd)
    pr.dump_stats(path)
    return path


def _merge_stats(paths: list[str]) -> pstats.Stats:
    st = pstats.Stats()
    for p in paths:
        st.add(p)
    return st


def sample(smiles: list[str], workers: int) -> tuple[pstats.Stats, int]:
    """Run molecules across a process pool and merge their profiles.

    Small inputs profile in-process (no spawn overhead); large inputs fan out
    across a process pool.
    """
    if len(smiles) < 300 or workers <= 1:
        _init_worker()
        path = _profile_chunk(smiles, progress=_report)
        try:
            return _merge_stats([path]), len(smiles)
        finally:
            try:
                os.unlink(path)
            except OSError:
                pass

    n_workers = min(max(2, workers), max(2, len(smiles) // 150))
    chunk_size = max(1, len(smiles) // (n_workers * 4))
    chunks = [smiles[i:i + chunk_size] for i in range(0, len(smiles), chunk_size)]
    paths: list[str] = []
    done_mols = 0
    try:
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker) as ex:
            fut_to_chunk = {ex.submit(_profile_chunk, ch): ch for ch in chunks}
            for fut in as_completed(fut_to_chunk):
                paths.append(fut.result())
                done_mols += len(fut_to_chunk[fut])
                _report(done_mols, len(smiles))
        return _merge_stats(paths), len(smiles)
    finally:
        for p in paths:
            try:
                os.unlink(p)
            except OSError:
                pass


def _build_graph(
    st: pstats.Stats, module: str, aggregate_external: bool
) -> tuple[list[dict], list[dict], float]:
    """Build nodes/edges from merged pstats (ported from the inline script)."""
    total_s = st.total_tt or 1e-9
    BS = chr(92)

    keep = {}
    for key, val in st.stats.items():
        fname, lineno, func = key
        if module not in fname.replace(BS, "/"):
            continue
        if func == "<module>":
            continue
        keep[key] = val

    maxct = max((v[3] for v in keep.values()), default=1e-9)
    id_of = {key: i for i, key in enumerate(keep)}

    nodes = []
    for key, (cc, nc, tt, ct, callers) in keep.items():
        fname, lineno, func = key
        rel = fname.replace(BS, "/").split(module + "/")[-1]
        m = re.match("layer([0-9])/", rel)
        layer = int(m.group(1)) if m else (6 if rel.startswith("tools/") else None)
        nodes.append(
            {
                "id": id_of[key],
                "label": func,
                "module": rel,
                "layer": layer,
                "line": lineno,
                "file": fname,
                "ncalls": cc,
                "self_s": round(tt, 6),
                "cum_s": round(ct, 6),
                "cum_pct": round(100.0 * ct / total_s, 3),
                "self_pct": round(100.0 * tt / total_s, 3),
                "w": round(ct / maxct, 4),
            }
        )

    edge_rows = []
    for key, (cc, nc, tt, ct, callers) in keep.items():
        for ck, (ccc, ncc, ttc, ctc) in callers.items():
            if ck in keep and ctc > 0:
                edge_rows.append((key, ck, ccc, ctc))
    max_ect = max((r[3] for r in edge_rows), default=1e-9)

    edges = []
    for (callee_key, caller_key, ccc, ctc) in edge_rows:
        edges.append(
            {
                "from": id_of[caller_key],
                "to": id_of[callee_key],
                "calls": ccc,
                "cum_s": round(ctc, 6),
                "cum_pct": round(100.0 * ctc / total_s, 3),
                "w": round(ctc / max_ect, 4),
            }
        )

    if aggregate_external:
        ext_ns = defaultdict(
            lambda: {"calls": 0, "cum_s": 0.0,
                     "callers": defaultdict(lambda: {"calls": 0, "cum_s": 0.0})}
        )
        for key, (cc, nc, tt, ct, callers) in st.stats.items():
            fname, lineno, func = key
            if module in fname.replace(BS, "/"):
                continue
            fn = fname.replace(BS, "/")
            parts = fn.split("/")
            ns = parts[-2] if len(parts) >= 2 else fn
            for i, seg in enumerate(parts):
                if seg in ("site-packages", "Lib") and i + 1 < len(parts):
                    ns = parts[i + 1]
                    break
            for ck, (ccc, ncc, ttc, ctc) in callers.items():
                if ck in keep and ctc > 0:
                    d = ext_ns[ns]
                    d["calls"] += ccc
                    d["cum_s"] += ctc
                    d["callers"][id_of[ck]]["calls"] += ccc
                    d["callers"][id_of[ck]]["cum_s"] += ctc
        nid = len(nodes)
        for ns, d in sorted(ext_ns.items()):
            if d["cum_s"] <= 0:
                continue
            nodes.append(
                {
                    "id": nid,
                    "label": ns,
                    "module": ns,
                    "layer": None,
                    "line": None,
                    "file": ns,
                    "ncalls": d["calls"],
                    "self_s": round(d["cum_s"], 6),
                    "cum_s": round(d["cum_s"], 6),
                    "cum_pct": round(100.0 * d["cum_s"] / total_s, 3),
                    "self_pct": round(100.0 * d["cum_s"] / total_s, 3),
                    "w": round(d["cum_s"] / maxct, 4),
                }
            )
            for cid, cd in d["callers"].items():
                edges.append(
                    {
                        "from": cid,
                        "to": nid,
                        "calls": cd["calls"],
                        "cum_s": round(cd["cum_s"], 6),
                        "cum_pct": round(100.0 * cd["cum_s"] / total_s, 3),
                        "w": round(cd["cum_s"] / max_ect, 4),
                    }
                )
            nid += 1

    return nodes, edges, total_s


def main() -> None:
    ap = argparse.ArgumentParser(description="Multi-process cProfile sampler")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--module", default="namepredict")
    ap.add_argument("--agg", action="store_true")
    ap.add_argument("--data", default=str(ROOT / "data" / "merged_benchmark.json"))
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()

    rows = json.load(open(args.data, encoding="utf-8"))
    smiles = [str(r["smiles"]) for r in rows[: args.n]]
    workers = args.workers or max(1, (os.cpu_count() or 1) - 1)

    t0 = time.perf_counter()
    st, n_actual = sample(smiles, workers)
    wall_s = time.perf_counter() - t0
    nodes, edges, total_s = _build_graph(st, args.module, args.agg)

    print("__CALLGRAPH_JSON__")
    print(
        json.dumps(
            {
                "nodes": nodes,
                "edges": edges,
                "total_s": total_s,  # CPU 累计（占比分母）
                "wall_s": round(wall_s, 3),  # 墙钟采样耗时
                "n_actual": n_actual,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
