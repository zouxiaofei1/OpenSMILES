"""Parallel bilingual benchmark scorer (faster; does not replace sequential script).

Usage:
  python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json
  python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --workers 8
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_str = str(_SRC)
if _src_str not in sys.path:
    sys.path.insert(0, _src_str)

from benchmarks.benchmark import (  # noqa: E402
    _SNAPSHOT_PATH,
    _empty_bucket,
    _fail_entry,
    _finalize,
    _handle_report,
    _is_fail,
    _load_rows,
    _print_summary,
    _result_entry,
    _tally,
    score_record,
)

# Per-process namer (set by worker initializer; not shared across processes).
_WORKER_NAMER: Any = None


def _init_worker() -> None:
    """Create one SMILESNNamer per worker process (avoids pickling namer)."""
    global _WORKER_NAMER
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _score_row(row: dict[str, Any]) -> tuple[dict[str, Any], str, str, dict[str, Any]]:
    """Score one gold row in a worker. Returns (score, pred_en, pred_zh, row)."""
    result = _WORKER_NAMER.name(str(row.get("smiles") or ""))
    pred_en, pred_zh = result.en or "", result.zh or ""
    score = score_record(pred_en, pred_zh, row)
    return score, pred_en, pred_zh, row


def _default_workers() -> int:
    n = os.cpu_count() or 1
    # Leave one core free on multi-core machines; at least 1.
    return max(1, n - 1) if n > 1 else 1


def _chunksize(n_rows: int, workers: int) -> int:
    """Balance IPC overhead vs load distribution for ~4k rows."""
    if n_rows <= 0 or workers <= 0:
        return 1
    # Aim for ~4–8 task batches per worker.
    return max(1, n_rows // (workers * 6))


def _aggregate(
    results: list[tuple[dict[str, Any], str, str, dict[str, Any]]],
) -> dict[str, Any]:
    total = _empty_bucket()
    by_source: dict[str, dict[str, int]] = defaultdict(_empty_bucket)
    by_tier: dict[Any, dict[str, int]] = defaultdict(_empty_bucket)
    fails: list[dict] = []
    entries: list[dict] = []
    for score, pred_en, pred_zh, row in results:
        _tally(total, score)
        _tally(by_source[str(row.get("source") or "unknown")], score)
        _tally(by_tier[row.get("tier", 0)], score)
        entry = _result_entry(row, score, pred_en, pred_zh)
        entries.append(entry)
        if _is_fail(score):
            fails.append(entry)
    return _finalize(total, by_source, by_tier, fails, entries)


def run_benchmark_parallel(
    data_path: str | Path,
    limit: int | None = None,
    workers: int | None = None,
) -> dict[str, Any]:
    """Run namer on gold JSON with a process pool; same report shape as sequential."""
    rows = _load_rows(Path(data_path), limit)
    n_workers = max(1, workers if workers is not None else _default_workers())
    if n_workers == 1 or len(rows) <= 1:
        # Sequential path reuses one process (no pool overhead).
        _init_worker()
        return _aggregate([_score_row(r) for r in rows])

    chunksize = _chunksize(len(rows), n_workers)
    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker,
    ) as pool:
        results = list(pool.map(_score_row, rows, chunksize=chunksize))
    return _aggregate(results)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Parallel bilingual SMILES namer benchmark (faster alternative)",
    )
    p.add_argument("--data", type=Path, required=True, help="Path to merged_benchmark.json")
    p.add_argument("--limit", type=int, default=None, help="Optional row limit")
    p.add_argument("--json", action="store_true", help="Print full report as JSON")
    p.add_argument(
        "--workers",
        type=int,
        default=None,
        help=f"Process pool size (default: cpu_count-1 = {_default_workers()})",
    )
    p.add_argument(
        "--time",
        action="store_true",
        help="Print wall-clock seconds on stderr",
    )
    p.add_argument(
        "--snapshot",
        type=Path,
        default=_SNAPSHOT_PATH,
        help="Path to last-run snapshot for diff (default: benchmarks/.last_run.json)",
    )
    p.add_argument(
        "--no-snapshot",
        action="store_true",
        help="Skip load/save of last-run snapshot (no diff)",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    t0 = time.perf_counter()
    try:
        report = run_benchmark_parallel(
            args.data,
            limit=args.limit,
            workers=args.workers,
        )
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
    elapsed = time.perf_counter() - t0
    n_workers = args.workers if args.workers is not None else _default_workers()
    n_rows = int(report.get("n_dual") or report.get("n_en") or 0)
    avg_ms = (elapsed * 1000.0 / n_rows) if n_rows else 0.0
    if args.no_snapshot:
        if args.json:
            out = {k: v for k, v in report.items() if k != "results"}
            out["elapsed_sec"] = round(elapsed, 3)
            out["workers"] = n_workers
            out["avg_ms_per_row"] = round(avg_ms, 3)
            print(json.dumps(out, ensure_ascii=False))
        else:
            _print_summary(report)
        if args.time:
            print(
                f"workers={n_workers} elapsed={elapsed:.2f}s "
                f"avg={avg_ms:.2f}ms/row (n={n_rows})",
                file=sys.stderr,
            )
        return

    # Attach timing fields before handle so --json consumers still get them.
    report = dict(report)
    report["elapsed_sec"] = round(elapsed, 3)
    report["workers"] = n_workers
    report["avg_ms_per_row"] = round(avg_ms, 3)
    _handle_report(report, as_json=args.json, snapshot_path=args.snapshot)
    if args.time and not args.json:
        print(
            f"workers={n_workers} elapsed={elapsed:.2f}s "
            f"avg={avg_ms:.2f}ms/row (n={n_rows})",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
