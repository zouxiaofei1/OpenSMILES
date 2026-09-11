"""Parallel generator for the benchmark *preview* cache (per-row predictions).

Spawned as `python -m benchmarks.benchmark_preview_parallel` by
server/backend/routes_benchmark.py (live + historical-commit branches), replacing the
old single-process `for`-loop inline scripts. Naming each gold row is
CPU-heavy, so rows run on a process pool (one SMILESNNamer per worker, cache
cleared per row for determinism, mirroring benchmark_parallel.py).

Progress/ordering contract with the server — keep these stable:
  * The on-disk cache is always a *continuous prefix* of source rows, so the
    front end can render partial rows before the run finishes.
  * Every ~200 rows the cache is rewritten and stdout prints
    `PROGRESS <done>/<total> <rate>r/s` (the server's status endpoint reads
    cache length; the history job registry parses this line's counter).

Usage:
  python -m benchmarks.benchmark_preview_parallel --data <source.json> --cache <out.json>
      [--sig]          wrap output as {"_data_sig": "...", "rows": [...]} (history caches)
      [--src <dir>]    import namepredict from <dir> instead of ROOT/src (history worktree)
      [--workers N]    pool size; 0 = cpu_count - 1
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
_FLUSH_EVERY = 200  # cache rewrite cadence, matches the old inline generator

# One namer per worker process; never pickled/shared across processes.
_WORKER_NAMER: Any = None


def _default_workers() -> int:
    n = os.cpu_count() or 1
    return max(1, n - 1) if n > 1 else 1


def _init_worker() -> None:
    """Create one SMILESNNamer per worker; silence stdout/stderr + RDKit logs."""
    global _WORKER_NAMER
    # The main process owns the pipe that the server reads for PROGRESS lines;
    # worker chatter (incl. leftover debug prints in the naming pipeline) must
    # not interleave with it or flood the pipe buffer.
    devnull = open(os.devnull, "w", encoding="utf-8")
    sys.stdout = devnull
    sys.stderr = devnull
    try:
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.ERROR)
        # RDKit C++ logs (e.g. kekulize warnings) bypass the Python-level
        # sys.stderr redirect above — DisableLog silences them at the source.
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _derived(smi: str, pe: str, pz: str, ge: str, gz: str) -> dict[str, Any]:
    """稠环 / 复杂度 / 中英相似度; 与 server 的旧缓存补算共用同一实现。"""
    from benchmarks.preview_metrics import row_derived

    return row_derived(smi, pe, pz, ge, gz)


def _score(pe: str, pz: str, row: dict[str, Any]) -> dict[str, Any]:
    """Score one prediction vs the gold row (mirrors old inline score_pred).

    显式 eval_en/eval_zh 标记优先; eval_zh=False 的行不考核中文(ok 只取决于
    英文)。无 eval 字段的源 (chebi20_test_1k.json 等) 回退旧行为: gold 非空即考核。
    """
    from namepredict.constants import nospace, normalize_en, normalize_zh

    ge = str(row.get("english_name") or "")
    gz = str(row.get("chinese_name") or "")
    ee = row.get("eval_en")
    ez = row.get("eval_zh")
    use_en = bool(ge) if ee is None else bool(ee)
    use_zh = bool(gz) if ez is None else bool(ez)
    # 空格不敏感：与 benchmark.py 判分一致，分词空格差异不判分
    en_ok = None if not use_en else nospace(normalize_en(pe)) == nospace(normalize_en(ge))
    zh_ok = None if not use_zh else nospace(normalize_zh(pz)) == nospace(normalize_zh(gz))
    if en_ok is not None and zh_ok is not None:
        dual = en_ok and zh_ok
    elif en_ok is not None:
        dual = en_ok
    elif zh_ok is not None:
        dual = zh_ok
    else:
        dual = False
    return {"en_ok": en_ok, "zh_ok": zh_ok, "ok": dual, "ret": bool(pe or pz)}


def _process_row(row: dict[str, Any]) -> dict[str, Any]:
    """Name one gold row in a worker and return its preview payload row.

    Clear the worker namer cache before each row: the naming pipeline's cache
    is stateful, and completion order under a pool is unpredictable — clearing
    restores the same per-row determinism benchmark_parallel relies on.
    """
    smi = str(row.get("smiles") or "")
    ge = str(row.get("english_name") or "")
    gz = str(row.get("chinese_name") or "")
    try:
        _WORKER_NAMER.cache.clear()
        result = _WORKER_NAMER.name(smi)
        pe, pz = result.en or "", result.zh or ""
    except Exception:
        pe, pz = "", ""
    return {
        "s": smi,
        "en": pe,
        "zh": pz,
        "ge": ge,
        "gz": gz,
        **_derived(smi, pe, pz, ge, gz),
        **_score(pe, pz, row),
    }


def _empty_payload(row: dict[str, Any]) -> dict[str, Any]:
    smi = str(row.get("smiles") or "")
    ge = str(row.get("english_name") or "")
    gz = str(row.get("chinese_name") or "")
    return {
        "s": smi,
        "en": "",
        "zh": "",
        "ge": ge,
        "gz": gz,
        **_derived(smi, "", "", ge, gz),
        **_score("", "", row),
    }


def _atomic_write(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--data", type=Path, required=True, help="Source benchmark rows (.json list)")
    p.add_argument("--cache", type=Path, required=True, help="Output cache file to write")
    p.add_argument("--sig", action="store_true",
                   help="Wrap output as {\"_data_sig\": ..., \"rows\": [...]} (history caches)")
    p.add_argument("--src", type=Path, default=None,
                   help="Import namepredict from this dir (default: ROOT/src)")
    p.add_argument("--workers", type=int, default=0,
                   help=f"Process pool size (default: cpu_count-1 = {_default_workers()})")
    return p


def _main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    source: Path = args.data
    cache: Path = args.cache
    src = args.src if args.src is not None else ROOT / "src"

    # Must precede any namepredict import *and* any pool spawn: worker children
    # inherit this sys.path (spawn), so imports there resolve to the same source.
    src_str = str(src)
    if src_str not in sys.path:
        sys.path.insert(0, src_str)
    try:
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.ERROR)
        # RDKit C++ logs (e.g. kekulize warnings) bypass the Python-level
        # sys.stderr redirect above — DisableLog silences them at the source.
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass

    rows = json.loads(source.read_text(encoding="utf-8"))
    total = len(rows)
    workers = args.workers if args.workers > 0 else _default_workers()
    sig = None
    if args.sig:
        try:
            st = source.stat()
            sig = f"{st.st_size}:{st.st_mtime_ns}"
        except OSError:
            sig = "missing"

    t0 = time.perf_counter()
    if total == 0:
        _atomic_write(cache, {"_data_sig": sig, "rows": []} if args.sig else [])
        print(f"DONE 0 rows 0.0s dual_ok=0", flush=True)
        return 0

    # ordered[prefix:] holds finished rows by source index; only the contiguous
    # prefix is ever flushed, so the cache never contains holes.
    ordered: list[dict[str, Any] | None] = [None] * total
    prefix = 0
    flushed = 0
    dual_ok = 0

    def write_prefix() -> None:
        nonlocal flushed
        payload = [r for r in ordered[:prefix] if r is not None]
        _atomic_write(cache, {"_data_sig": sig, "rows": payload} if args.sig else payload)
        flushed = prefix
        elapsed = time.perf_counter() - t0
        rate = (prefix / elapsed) if elapsed > 0 else 0.0
        print(f"PROGRESS {prefix}/{total} {rate:.1f}r/s", flush=True)

    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        futures = {pool.submit(_process_row, row): i for i, row in enumerate(rows)}
        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                ordered[idx] = fut.result()
            except Exception:
                # A worker died (e.g. native crash): score the row as empty rather
                # than dropping it, so every source index still yields a payload.
                ordered[idx] = _empty_payload(rows[idx])
            if ordered[idx] is not None and ordered[idx].get("ok"):
                dual_ok += 1
            while prefix < total and ordered[prefix] is not None:
                prefix += 1
            if prefix - flushed >= _FLUSH_EVERY or prefix == total:
                write_prefix()

    elapsed = time.perf_counter() - t0
    print(f"DONE {total} rows {elapsed:.1f}s dual_ok={dual_ok}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(_main())
