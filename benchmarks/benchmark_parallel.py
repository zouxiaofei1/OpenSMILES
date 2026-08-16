"""Parallel bilingual benchmark scorer (faster; does not replace sequential script).

Usage:
  python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json
  python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --workers 8
"""

from __future__ import annotations
try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except Exception:
    pass

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_str = str(_SRC)
if _src_str not in sys.path:
    sys.path.insert(0, _src_str)

from benchmarks.benchmark import (  # noqa: E402
    _SNAPSHOT_PATH,
    _empty_bucket,
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
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _score_row(row: dict[str, Any]) -> tuple[dict[str, Any], str, str, dict[str, Any]]:
    """Score one gold row in a worker. Returns (score, pred_en, pred_zh, row).

    每行命名前清空 worker 共享 cache：命名管线依赖 cache 状态（同一 SMILES 的
    递归子结构命中/未命中路径不同，结果可能不同），行分配/处理顺序在并行下不可
    预测，共享 cache 会让同名结果跨行漂移。每行从一致初始状态命名保证确定性。
    """
    try:
        _WORKER_NAMER.cache.clear()
        result = _WORKER_NAMER.name(str(row.get("smiles") or ""))
        pred_en, pred_zh = result.en or "", result.zh or ""
    except Exception:
        pred_en, pred_zh = "", ""
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


_PROGRESS_EVERY = 100


def _progress_bar(done: int, total: int, width: int = 28) -> str:
    if total <= 0:
        return f"[{'?' * width}] 0/0"
    frac = min(1.0, done / total)
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {done}/{total} ({100.0 * frac:.1f}%)"


def _fmt_acc(ok: int, n: int) -> str:
    pct = 0.0 if n == 0 else 100.0 * ok / n
    return f"{ok}/{n} ({pct:.1f}%)"


def _print_progress(done: int, total: int, t0: float, bucket: dict[str, int]) -> None:
    """Print bar + running EN/ZH/dual accuracy so far (preview)."""
    elapsed = time.perf_counter() - t0
    rate = (done / elapsed) if elapsed > 0 else 0.0
    eta = ((total - done) / rate) if rate > 0 else 0.0
    en = _fmt_acc(bucket["ok_en"], bucket["n_en"])
    zh = _fmt_acc(bucket["ok_zh"], bucket["n_zh"])
    dual = _fmt_acc(bucket["ok_dual"], bucket["n_dual"])
    print(
        f"progress {_progress_bar(done, total)}  "
        f"{rate:.1f} row/s  elapsed={elapsed:.1f}s  eta={eta:.1f}s\n"
        f"  so far  EN={en}  ZH={zh}  dual={dual}",
        flush=True,
    )


def _collect_with_progress(iterable, total: int, t0: float | None = None):
    """Drain scored rows; every 100 print bar + running EN/ZH/dual accuracy."""
    t0 = time.perf_counter() if t0 is None else t0
    results: list[tuple[dict[str, Any], str, str, dict[str, Any]]] = []
    bucket = _empty_bucket()
    for i, item in enumerate(iterable, 1):
        results.append(item)
        _tally(bucket, item[0])
        if i % _PROGRESS_EVERY == 0 or i == total:
            _print_progress(i, total, t0, bucket)
    return results


def _print_timeouts(indexes: list[int], timeout: float) -> None:
    if not indexes:
        return
    joined = ",".join(str(index) for index in indexes)
    print(f"timeout index={joined} after={timeout:g}s", file=sys.stderr, flush=True)


def _pool_results_with_progress(
    rows: list[dict[str, Any]],
    workers: int,
    timeout: float,
    t0: float,
) -> list[tuple[dict[str, Any], str, str, dict[str, Any]]]:
    """Run rows on persistent process pool with per-row timeout; report progress."""
    total = len(rows)
    results: list[tuple[dict[str, Any], str, str, dict[str, Any]]] = []
    bucket = _empty_bucket()
    timeout_indexes: list[int] = []

    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        futures = {pool.submit(_score_row, row): i for i, row in enumerate(rows)}
        for i, fut in enumerate(as_completed(futures), 1):
            idx = futures[fut]
            try:
                item = fut.result(timeout=timeout)
            except Exception:
                row = rows[idx]
                score = score_record("", "", row)
                item = (score, "", "", row)
                timeout_indexes.append(idx)
            results.append(item)
            _tally(bucket, item[0])
            if i % _PROGRESS_EVERY == 0 or i == total:
                _print_progress(i, total, t0, bucket)
    _print_timeouts(timeout_indexes, timeout)
    return results


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
    timeout: float = 1.0,
) -> dict[str, Any]:
    """Run namer on gold JSON with a process pool; same report shape as sequential."""
    rows = _load_rows(Path(data_path), limit)
    n_rows = len(rows)
    n_workers = max(1, workers if workers is not None else _default_workers())
    t0 = time.perf_counter()
    print(f"benchmark start: n={n_rows} workers={n_workers}", flush=True)
    if n_rows == 0:
        return _aggregate([])

    results = _pool_results_with_progress(rows, n_workers, timeout, t0)
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
        "--timeout",
        type=float,
        default=1.0,
        help="Per-row hard timeout in seconds (default: 1.0)",
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
            timeout=args.timeout,
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
