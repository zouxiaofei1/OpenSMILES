"""OPSIN 双向翻译基准：SMILES → namer → OPSIN 回译 → 结构一致才算通过。

与 benchmark.py 的字符串打分不同，本基准的 en 维只看回译结构是否与原 SMILES
一致（MATCH 通过，STEREO_ONLY / DIFF / OPSIN_FAIL / NAMER_FAIL 不通过）。
数据格式、报告、快照、diff 与 CLI 全部复用 benchmark.py / benchmark_parallel.py
那一套：同 schema 的 opsin_benchmark.json，同 en/zh/dual 报告与 IMPROVE/REGRESS 差分。

用法:
  python -m benchmarks.opsin_benchmark --data benchmarks/opsin_benchmark.json
  python -m benchmarks.opsin_benchmark --data benchmarks/opsin_benchmark.json --time --timeout 1
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
import tempfile
import time
import warnings
from collections import Counter, defaultdict
from concurrent.futures import TimeoutError as FutureTimeout
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from pebble import ProcessPool

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from benchmarks.benchmark import (  # noqa: E402
    _empty_bucket,
    _finalize,
    _handle_report,
    _load_rows,
    _print_summary,
    _result_entry,
    _tally,
)

# Last-run snapshot for diff vs previous pass/verdict results.
_SNAPSHOT_PATH = Path(__file__).resolve().parent / ".opsin_last_run.json"
OPSIN_BATCH = 2000
PASS_VERDICTS = {"MATCH"}

# OPSIN 只吃 ASCII：把上下标/希腊字母等记法折平。
_TRANS = str.maketrans({
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6",
    "⁷": "7", "⁸": "8", "⁹": "9", "⁺": "+", "⁻": "-", "⁼": "=",
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6",
    "₇": "7", "₈": "8", "₉": "9",
    "λ": "lambda", "Λ": "lambda", "μ": "mu", "α": "alpha", "β": "beta",
    "γ": "gamma", "δ": "delta", "ω": "omega", "κ": "kappa",
    "′": "'", "″": '"', "·": ".", "×": "x", "–": "-", "—": "-",
    "’": "'", "‘": "'", "\xa0": " ",
})

_WORKER_NAMER: Any = None


def _norm_ascii(name: str) -> str:
    """非 ASCII 记号转 ASCII，去掉控制字符。"""
    s = (name or "").translate(_TRANS)
    return "".join(c for c in s if ord(c) >= 32).strip()


def _row_key(row: dict[str, Any]) -> str:
    """取行标识：id 优先，回落 smiles（与 benchmark.py 一致）。"""
    rid = row.get("id")
    if rid is not None and str(rid) != "":
        return f"id:{rid}"
    return f"smiles:{row.get('smiles') or ''}"


def _init_worker() -> None:
    """worker 初始化：静音、建 namer（每个进程一个）。"""
    global _WORKER_NAMER
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    from opensmiles.namer import SMILESNNamer
    _WORKER_NAMER = SMILESNNamer()


def _name_row(row: dict[str, Any]) -> tuple[str, bool, str]:
    """命名一条，返回 (key, namer_ok, name_en)。"""
    key = _row_key(row)
    _WORKER_NAMER.cache.clear()
    try:
        result = _WORKER_NAMER.name(str(row.get("smiles") or ""))
        return key, bool(result.success), (result.en or "")
    except Exception:
        return key, False, ""


def _opsin_task(args: tuple[int, list[str]]) -> tuple[int, list[str]]:
    """独立 JVM 批量回译一批名字，返回 (下标偏移, SMILES 列表)。"""
    off, names = args
    os.chdir(tempfile.mkdtemp(prefix="opsin_"))
    warnings.simplefilter("ignore")
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    from py2opsin import py2opsin
    try:
        return off, py2opsin([_norm_ascii(n) for n in names], output_format="SMILES",
                             allow_acid=True, allow_radicals=True, allow_bad_stereo=True)
    except Exception:
        return off, [""] * len(names)


def _verdict(smiles: str, namer_ok: bool, name_en: str, back: str) -> str:
    """按回译结构判定：MATCH / STEREO_ONLY / DIFF / OPSIN_FAIL / NAMER_FAIL / PARSE_FAIL。"""
    from rdkit import Chem
    if not namer_ok or not name_en:
        return "NAMER_FAIL"
    if not back:
        return "OPSIN_FAIL"
    a, b = Chem.MolFromSmiles(smiles), Chem.MolFromSmiles(back)
    if a is None or b is None:
        return "PARSE_FAIL"
    if Chem.MolToSmiles(a) == Chem.MolToSmiles(b):
        return "MATCH"
    if Chem.MolToSmiles(a, isomericSmiles=False) == Chem.MolToSmiles(b, isomericSmiles=False):
        return "STEREO_ONLY"
    return "DIFF"


def _default_workers() -> int:
    """默认 worker 数：核数减一。"""
    n = os.cpu_count() or 1
    return max(1, n - 1) if n > 1 else 1


def _name_all(rows: list[dict[str, Any]], workers: int, timeout: float) -> dict[str, tuple]:
    """多进程命名，返回 row_key -> (namer_ok, name_en, err)。"""
    t0 = time.perf_counter()
    out: dict[str, tuple] = {}
    pool = ProcessPool(max_workers=workers, initializer=_init_worker)
    try:
        it = iter(pool.map(_name_row, rows, timeout=timeout, chunksize=1).result())
        k = 0
        while True:
            try:
                key, ok, en = next(it)
                out[key] = (ok, en, "")
            except StopIteration:
                break
            except FutureTimeout:
                out[_row_key(rows[k])] = (False, "", "Timeout")
            except Exception as exc:
                out[_row_key(rows[k])] = (False, "", type(exc).__name__)
            k += 1
            if k % 200 == 0 or k == len(rows):
                el = time.perf_counter() - t0
                print(f"  命名 {k}/{len(rows)} {k / max(el, 1e-9):.0f}条/s",
                      file=sys.stderr, flush=True)
    finally:
        pool.stop()
        pool.join()
    return out


def _opsin_all(rows: list[dict[str, Any]], named: dict[str, tuple], workers: int) -> dict[str, str]:
    """批量回译，返回 row_key -> back_smiles。"""
    todo = [(i, named[_row_key(r)][1]) for i, r in enumerate(rows) if named[_row_key(r)][1]]
    if not todo:
        return {}
    tasks = [(i, [x[1] for x in todo[i:i + OPSIN_BATCH]])
             for i in range(0, len(todo), OPSIN_BATCH)]
    backs: dict[str, str] = {}
    with Pool(processes=workers, maxtasksperchild=1) as pool:
        for off, smis in pool.imap_unordered(_opsin_task, tasks):
            for j, s in enumerate(smis):
                backs[_row_key(rows[todo[off + j][0]])] = s or ""
    return backs


def run_benchmark(
    data_path: str | Path,
    limit: int | None = None,
    workers: int | None = None,
    timeout: float = 1.0,
) -> dict[str, Any]:
    """跑一遍双向翻译基准，返回与 benchmark.py 同形的报告。"""
    rows = _load_rows(Path(data_path), limit)
    n_workers = max(1, workers if workers is not None else _default_workers())
    print(f"opsin benchmark start: n={len(rows)} workers={n_workers}",
          file=sys.stderr, flush=True)

    named = _name_all(rows, n_workers, timeout)
    backs = _opsin_all(rows, named, n_workers)

    total = _empty_bucket()
    by_source: dict[str, dict[str, int]] = defaultdict(_empty_bucket)
    by_tier: dict[Any, dict[str, int]] = defaultdict(_empty_bucket)
    fails: list[dict] = []
    results: list[dict] = []
    verdicts: Counter = Counter()
    for row in rows:
        key = _row_key(row)
        ok, name_en, _err = named.get(key, (False, "", "Missing"))
        verdict = _verdict(str(row.get("smiles") or ""), ok, name_en, backs.get(key, ""))
        verdicts[verdict] += 1
        # en 维 = 回译结构是否一致；无 gold 中文名，zh 维不考核（eval_zh=false）。
        en_ok = (verdict in PASS_VERDICTS) if row.get("eval_en", True) else None
        score = {"en_ok": en_ok, "zh_ok": None, "dual_ok": bool(en_ok)}
        _tally(total, score)
        _tally(by_source[str(row.get("source") or "unknown")], score)
        _tally(by_tier[row.get("tier", 0)], score)
        entry = _result_entry(row, score, name_en, "")
        results.append(entry)
        if not en_ok and en_ok is not None:
            fails.append(entry)

    report = _finalize(total, by_source, by_tier, fails, results)
    report["verdicts"] = dict(verdicts.most_common())
    return report


def _build_parser() -> argparse.ArgumentParser:
    """命令行参数，对齐 benchmark_parallel.py。"""
    p = argparse.ArgumentParser(description="OPSIN round-trip benchmark (namer → OPSIN → structure)")
    p.add_argument("--data", type=Path, required=True, help="Path to benchmarks/opsin_benchmark.json")
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
        help="Per-row naming timeout in seconds (default: 1.0)",
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
        help="Path to last-run snapshot for diff (default: benchmarks/.opsin_last_run.json)",
    )
    p.add_argument(
        "--no-snapshot",
        action="store_true",
        help="Skip load/save of last-run snapshot (no diff)",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    """解析参数并跑基准。"""
    args = _build_parser().parse_args(argv)
    t0 = time.perf_counter()
    try:
        report = run_benchmark(
            args.data, limit=args.limit, workers=args.workers, timeout=args.timeout,
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
