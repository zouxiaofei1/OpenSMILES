"""Parallel stereo benchmark (faster; sequential in benchmarks/stereo_benchmark.py).

只比名称里的立体描述符 token（R/S/E/Z，含 fused 位次如 3aR），忽略整名其它差异。
与 benchmark_parallel 同构：worker 每行清缓存保证确定性、进度条、per-row timeout、
快照 diff（默认 benchmarks/.last_stereo.json）。

Usage:
  python -m benchmarks.stereo_benchmark_parallel --data benchmarks/merged_benchmark.json
  python -m benchmarks.stereo_benchmark_parallel --data benchmarks/merged_benchmark.json --workers 8
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_str = str(_SRC)
if _src_str not in sys.path:
    sys.path.insert(0, _src_str)

# 复用串行版的 token 抽取与快照/汇总/差异打印（口径一致）
from benchmarks import stereo_benchmark as _sb  # noqa: E402

# Per-process namer（worker 里建，不进 pickle）。
_WORKER_NAMER: Any = None


def _init_worker() -> None:
    """每个 worker 一个命名器并静音 src 调试 print；进度由主进程打印。"""
    global _WORKER_NAMER
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    from namepredict.namer import SMILESNNamer

    _WORKER_NAMER = SMILESNNamer()


def _name_row_worker(row: dict[str, Any]) -> dict[str, Any]:
    """在 worker 命名单行并返回立体比对字段（异常记 errored，不抛）。

    每行清 worker 共享 cache：行分配顺序不可预测，跨行共享片段缓存会让结果漂移；
    从一致初始状态命名保证并行确定性。
    """
    ref_t = _sb.stereo_tokens(row.get("english_name") or "")
    out = {"ref_tokens": sorted(ref_t), "our_tokens": [], "pred_en": "",
           "errored": False, "stereo_equal": False, "count_equal": False, "flip_free": False}
    if not ref_t:
        return out
    try:
        _WORKER_NAMER.cache.clear()
        rr = _WORKER_NAMER.name(str(row.get("smiles") or ""))
    except Exception:
        out["errored"] = True
        return out
    our_t = _sb.stereo_tokens(rr.en if rr.success else "")
    out["our_tokens"] = sorted(our_t)
    out["pred_en"] = rr.en if rr.success else ""
    if not (rr.success and rr.en and our_t):
        return out
    equal = Counter(our_t) == Counter(ref_t)
    out["stereo_equal"] = equal
    if len(our_t) == len(ref_t):
        out["count_equal"] = True
        out["flip_free"] = equal
    return out


def _default_workers() -> int:
    n = os.cpu_count() or 1
    return max(1, n - 1) if n > 1 else 1


_PROGRESS_EVERY = 100


def _progress_bar(done: int, total: int, width: int = 28) -> str:
    if total <= 0:
        return f"[{'?' * width}] 0/0"
    frac = min(1.0, done / total)
    filled = int(width * frac)
    return f"[{'#' * filled}{'-' * (width - filled)}] {done}/{total} ({100.0 * frac:.1f}%)"


def _fmt_acc(ok: int, n: int) -> str:
    pct = 0.0 if n == 0 else 100.0 * ok / n
    return f"{ok}/{n} ({pct:.1f}%)"


def _print_progress(done: int, total: int, t0: float, s: dict[str, int]) -> None:
    """打印进度条 + 迄今立体计数。"""
    elapsed = time.perf_counter() - t0
    rate = (done / elapsed) if elapsed > 0 else 0.0
    eta = ((total - done) / rate) if rate > 0 else 0.0
    equal = _fmt_acc(s["stereo_equal"], s["our_stereo"])
    pure = s["flip_free"] / (s["count_equal"] or 1)
    print(
        f"progress {_progress_bar(done, total)}  {rate:.1f} row/s  "
        f"elapsed={elapsed:.1f}s  eta={eta:.1f}s\n"
        f"  so far  stereo_equal={equal}  pure_ok={pure:.1%}  "
        f"count_equal={s['count_equal']}",
        flush=True,
    )


def _run_pool(rows: list[dict[str, Any]], workers: int, timeout: float | None) -> dict[str, Any]:
    """在进程池上命名全部行，边收边统计；返回 {stats, items} 供快照 diff。"""
    total = len(rows)
    stats = {"rows": total, "named_ok": 0, "errored": 0, "our_stereo": 0,
             "stereo_equal": 0, "count_equal": 0, "flip_free": 0}
    items: dict[str, dict[str, Any]] = {}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        futures = {pool.submit(_name_row_worker, row): i for i, row in enumerate(rows)}
        for i, fut in enumerate(as_completed(futures), 1):
            idx = futures[fut]
            try:
                res = fut.result(timeout=timeout)
            except Exception:
                res = {"ref_tokens": _sb.stereo_tokens(rows[idx].get("english_name") or ""),
                       "our_tokens": [], "errored": True,
                       "stereo_equal": False, "count_equal": False, "flip_free": False}
            key = _sb._row_key(rows[idx], idx)
            items[key] = res
            if res["errored"]:
                stats["errored"] += 1
                continue
            if res["our_tokens"]:
                stats["our_stereo"] += 1
                stats["named_ok"] += 1
            if res["stereo_equal"]:
                stats["stereo_equal"] += 1
            if res["count_equal"]:
                stats["count_equal"] += 1
                if res["flip_free"]:
                    stats["flip_free"] += 1
            if i % _PROGRESS_EVERY == 0 or i == total:
                _print_progress(i, total, t0, stats)
    stats["pure_ratio"] = stats["flip_free"] / (stats["count_equal"] or 1)
    return {"stats": stats, "items": items}


def _select_rows(data_path: Path, *, min_len: int, sample: int, limit: int | None) -> list[dict]:
    """选立体行（参考名含立体 token）；min_len 过滤名长，sample>0 取最长 sample 条。"""
    rows = _sb._load_rows(data_path)
    rows = [
        r for r in rows
        if (r.get("english_name") or "") and _sb.stereo_tokens(r.get("english_name") or "")
        and len(r.get("english_name") or "") >= min_len
    ]
    rows.sort(key=lambda r: -len(r.get("english_name") or ""))
    if sample:
        rows = rows[:sample]
    if limit is not None:
        rows = rows[:limit]
    return rows


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Parallel stereo benchmark (立体 token 正确比例；较串行更快)",
    )
    p.add_argument("--data", type=Path, required=True, help="Path to benchmarks/merged_benchmark.json")
    p.add_argument("--limit", type=int, default=None, help="Optional row limit")
    p.add_argument("--min-len", type=int, default=0, help="参考英文名最小长度")
    p.add_argument("--sample", type=int, default=0, help="取名长最长 N 条（0=全部立体行）")
    p.add_argument("--json", action="store_true", help="Print full report as JSON")
    p.add_argument("--workers", type=int, default=None,
                   help=f"Process pool size (default: cpu_count-1 = {_default_workers()})")
    p.add_argument("--timeout", type=float, default=None,
                   help="Per-row hard timeout in seconds (default: none)")
    p.add_argument("--time", action="store_true", help="Print wall-clock seconds on stderr")
    p.add_argument("--snapshot", type=Path, default=_sb._SNAPSHOT_PATH)
    p.add_argument("--no-snapshot", action="store_true", help="Skip load/save of last snapshot")
    p.add_argument("--errors", type=Path, default=None,
                   help="把立体不匹配/命名报错的行导出为 JSON")
    return p


def main(argv: list[str] | None = None) -> None:
    """CLI：并行测量立体正确比例，可写/比快照（与 benchmark_parallel 同风格）。"""
    args = _build_parser().parse_args(argv)
    t0 = time.perf_counter()
    try:
        rows = _select_rows(args.data, min_len=args.min_len, sample=args.sample, limit=args.limit)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
    n_workers = args.workers if args.workers is not None else _default_workers()
    print(f"stereo benchmark start: n={len(rows)} workers={n_workers}", flush=True)
    report = _run_pool(rows, n_workers, args.timeout)
    elapsed = time.perf_counter() - t0
    avg_ms = (elapsed * 1000.0 / len(rows)) if rows else 0.0
    if args.errors is not None:
        _sb.dump_errors(rows, report, args.errors)

    if args.no_snapshot:
        _sb._print_summary(report, elapsed)
    else:
        prev = _sb._load_snapshot(args.snapshot)
        _sb._print_summary(report, elapsed)
        if prev is None:
            print("diff_vs_last: none (first run)")
        else:
            _sb._print_diffs(_sb._diff(prev, report))
        _sb._save_snapshot(report, args.snapshot)
    if args.time:
        print(f"workers={n_workers} elapsed={elapsed:.2f}s avg={avg_ms:.2f}ms/row (n={len(rows)})",
              file=sys.stderr)


if __name__ == "__main__":
    main()
