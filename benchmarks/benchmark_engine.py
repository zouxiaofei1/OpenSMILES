"""Run the gold bilingual benchmark against an arbitrary naming engine dir.

Reuses the scorer/report helpers from benchmarks.benchmark (same normalize rules,
same by_source/by_tier buckets).  Purpose-built so an engine under tools/ (e.g.
tools/namepredict-v3, whose package dir has a hyphen and a dict-style .name()
result) can be scored head-to-head with the active src namer without touching the
stable harnesses (benchmark.py / benchmark_parallel.py / their snapshots).

Usage:
  python -m benchmarks.benchmark_engine --engine v3 --metric loose  # default benchmarks/merged_benchmark.json
  python -m benchmarks.benchmark_engine --engine src --metric loose # reference: active src namer
  python -m benchmarks.benchmark_engine --engine v3 --limit 200 --workers 8 --metric strict
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from rdkit import RDLogger  # noqa: F401  (before disabling)

RDLogger.DisableLog("rdApp.*")

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_str = str(_SRC)
if _src_str not in sys.path:
    sys.path.insert(0, _src_str)

from benchmarks.benchmark import (  # noqa: E402
    _empty_bucket,
    _finalize,
    _is_fail,
    _load_rows,
    _result_entry,
    _tally,
    score_record,
)
from namepredict.tools.re import normalize_en, normalize_zh  # noqa: E402

# ── loose（约定不敏感）归一 ─────────────────────────────────────
# 目的：gold 是按现役 src 引擎的书写约定 curated 的（位次必带、取代基与母体是否加连字符、
# 是否加括号等），逐字严格比对会把"化学等价但书写不同"误判为错。loose 把括号字符、
# 连字符类分隔符与 ASCII 位次数字全部剥掉，只留字母/汉字词干，再判相等。
_BRACKET_RE = re.compile(r"[()\[\]（）【】]")  # 剥括号（保留内部文字）
_DROP_CHARS = set("-–—,，;；· .'’\"") | set("()[]（）【】")


def _loose_token(name: str) -> str:
    """去括号/连字符类分隔符/ASCII 位次数字，返回只含词干字母与汉字的串。"""
    s = name
    s = _BRACKET_RE.sub("", s)
    s = re.sub(r"[0-9]+", "", s)
    return "".join(ch for ch in s if ch not in _DROP_CHARS)


def _loose_en(name: str) -> str:
    """英文 loose 规范形：先 normalize_en（小写、统一连字符/括号），再剥分隔符与位次。"""
    return _loose_token(normalize_en(name))


def _loose_zh(name: str) -> str:
    """中文 loose 规范形：先 normalize_zh（统一括号种类），再剥分隔符与位次。"""
    return _loose_token(normalize_zh(name))


def _check_loose_en(pred: str, gold: str) -> bool:
    return _loose_en(pred) == _loose_en(gold)


def _check_loose_zh(pred: str, gold: str) -> bool:
    return _loose_zh(pred) == _loose_zh(gold)


def score_record_loose(pred_en: str, pred_zh: str, row: dict[str, Any]) -> dict[str, Any]:
    """Loose (括号/连字符/位次不敏感) field-wise scoring."""
    en_ok = _check_loose_en(pred_en, row.get("english_name") or "") if row.get("eval_en") else None
    zh_ok = _check_loose_zh(pred_zh, row.get("chinese_name") or "") if row.get("eval_zh") else None
    checks = [x for x in (en_ok, zh_ok) if x is not None]
    dual_ok = all(checks) if checks else False
    return {"en_ok": en_ok, "zh_ok": zh_ok, "dual_ok": dual_ok}

# Worker engine; instantiated per-process by _init_worker (not shared/pickled).
_WORKER_NAMER: Any = None

_PROGRESS_EVERY = 100


# ---------------------------------------------------------------- engine loading

def _load_legacy_package(dirname: str, alias: str) -> Any:
    """Import tools/<dirname> under a valid alias so relative imports resolve.

    目录名带连字符(如 namepredict-v3)不能直接 import：以 __init__.py 为 spec 入口、
    submodule_search_locations 指向引擎目录，注册到 sys.modules 后包内相对导入
    (.layer0 等) 按 alias 包名正常解析。加载一次后从 sys.modules 直取。
    """
    import importlib.util

    if alias in sys.modules:
        return sys.modules[alias]
    d = os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools", dirname)
    )
    spec = importlib.util.spec_from_file_location(
        alias, os.path.join(d, "__init__.py"), submodule_search_locations=[d]
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def _install_v2_shim() -> None:
    """v2 内部 cache/common_names.py `from namepredict.core.capitalization import ...`。

    旧版顶层 namepredict(带 core 子包)已不存在; 现役 src/namepredict 无 core。把
    namepredict.core.capitalization shim 挂到已加载的 namepredict 包上, 函数指向
    v2 自带实现, 使 v2 与现役 src 能在同进程共存。
    """
    import types

    pkg = _load_legacy_package("namepredict-v2", "namepredict_v2")
    cap_mod = sys.modules.get(f"{pkg.__name__}.core.capitalization")
    if cap_mod is None:
        from importlib import import_module as _imp

        cap_mod = _imp(f"{pkg.__name__}.core.capitalization")
    np = sys.modules.get("namepredict")
    if np is None:  # 未加载现役 src: 建一个空壳父包承载 shim
        np = types.ModuleType("namepredict")
        np.__path__ = []
        sys.modules["namepredict"] = np
    core = sys.modules.get("namepredict.core")
    if core is None:
        core = types.ModuleType("namepredict.core")
        core.__path__ = []
        sys.modules["namepredict.core"] = core
        np.core = core
    cap = sys.modules.get("namepredict.core.capitalization")
    if cap is None:
        cap = types.ModuleType("namepredict.core.capitalization")
        sys.modules["namepredict.core.capitalization"] = cap
        core.capitalization = cap
    cap.capitalize_chemical_name = cap_mod.capitalize_chemical_name


class _NoopCache:
    """v2 用全局单例缓存(实例无 cache 属性); 提供与 src/v3 一致的 clear() 空操作。"""

    def clear(self) -> None:
        return None


class _Engine:
    """Uniform engine facade: .name(smiles) -> {en, zh} object, .cache with clear()."""

    def __init__(self, namer: Any, style: str, cache: Any) -> None:
        self._namer = namer
        self._style = style
        self.cache = cache

    def name(self, smiles: str) -> SimpleNamespace:
        result = self._namer.name(smiles)
        if self._style == "dict":
            return SimpleNamespace(en=result.get("en") or "", zh=result.get("zh") or "")
        return SimpleNamespace(en=result.en or "", zh=result.zh or "")


def _make_engine(engine: str) -> _Engine:
    """Build the requested engine in this process."""
    if engine == "v2":
        _install_v2_shim()
        pkg = sys.modules["namepredict_v2"]
        namer = pkg.SMILESNNamerV2()
        return _Engine(namer, style="dict", cache=_NoopCache())
    if engine == "v3":
        pkg = _load_legacy_package("namepredict-v3", "namepredict_v3")
        namer = pkg.SMILESNNamerV2()
        return _Engine(namer, style="dict", cache=namer._cache)
    # engine == "src": active src namer (already on sys.path via _SRC)
    from namepredict.namer import SMILESNNamer

    namer = SMILESNNamer()
    return _Engine(namer, style="attr", cache=namer.cache)


# ---------------------------------------------------------------- worker & pool

def _init_worker() -> None:
    """One engine per worker process.  Engine chosen by env var so it survives Windows spawn."""
    global _WORKER_NAMER
    sys.stdout = open(os.devnull, "w", encoding="utf-8")
    engine = os.environ.get("NAMER_BENCH_ENGINE", "v3")
    _WORKER_NAMER = _make_engine(engine)


def _scorer_for(metric: str):
    return score_record if metric == "strict" else score_record_loose


def _score_row(row: dict[str, Any]) -> tuple[dict[str, Any], str, str, dict[str, Any]]:
    """Score one gold row in a worker.  Clears the shared cache first (see benchmark_parallel)."""
    metric = os.environ.get("NAMER_BENCH_METRIC", "strict")
    scorer = _scorer_for(metric)
    try:
        _WORKER_NAMER.cache.clear()
        result = _WORKER_NAMER.name(str(row.get("smiles") or ""))
        pred_en, pred_zh = result.en, result.zh
    except Exception:
        pred_en, pred_zh = "", ""
    score = scorer(pred_en, pred_zh, row)
    return score, pred_en, pred_zh, row


def _default_workers() -> int:
    n = os.cpu_count() or 1
    return max(1, n - 1) if n > 1 else 1


def _chunksize(n_rows: int, workers: int) -> int:
    if n_rows <= 0 or workers <= 0:
        return 1
    return max(1, n_rows // (workers * 6))


def _print_progress(done: int, total: int, t0: float, bucket: dict[str, int]) -> None:
    elapsed = time.perf_counter() - t0
    rate = (done / elapsed) if elapsed > 0 else 0.0
    print(
        f"  {done}/{total}  {rate:.1f} row/s  elapsed={elapsed:.0f}s  "
        f"en_ok={bucket['ok_en']}/{bucket['n_en']}  zh_ok={bucket['ok_zh']}/{bucket['n_zh']}",
        flush=True,
    )


def _run_pool(rows: list[dict[str, Any]], workers: int, timeout: float, t0: float):
    """Run rows on a persistent pool, return per-row score/pred items in submission order."""
    total = len(rows)
    out: list[tuple[dict, str, str, dict] | None] = [None] * total
    bucket = _empty_bucket()
    done = 0
    timeouts: list[int] = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        futures = {pool.submit(_score_row, row): i for i, row in enumerate(rows)}
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                item = fut.result(timeout=timeout)
            except Exception:
                row = rows[i]
                item = (score_record("", "", row), "", "", row)
                timeouts.append(i)
            out[i] = item
            done += 1
            _tally(bucket, item[0])
            if done % _PROGRESS_EVERY == 0 or done == total:
                _print_progress(done, total, t0, bucket)
    if timeouts:
        print(f"  timeout rows={len(timeouts)}", file=sys.stderr, flush=True)
    return [x for x in out if x is not None]


def _aggregate(results) -> dict[str, Any]:
    total = _empty_bucket()
    by_source: dict[str, dict[str, int]] = defaultdict(_empty_bucket)
    by_tier: dict[Any, dict[str, int]] = defaultdict(_empty_bucket)
    fails: list[dict] = []
    entries: list[dict] = []
    named_en = named_zh = 0
    for score, pred_en, pred_zh, row in results:
        _tally(total, score)
        _tally(by_source[str(row.get("source") or "unknown")], score)
        _tally(by_tier[row.get("tier", 0)], score)
        entry = _result_entry(row, score, pred_en, pred_zh)
        entries.append(entry)
        if _is_fail(score):
            fails.append(entry)
        if pred_en:
            named_en += 1
        if pred_zh:
            named_zh += 1
    report = _finalize(total, by_source, by_tier, fails, entries)
    report["named_en"] = named_en
    report["named_zh"] = named_zh
    return report


# ----------------------------------------------------------------------- report

def _fmt_pct(acc: float) -> str:
    return f"{100.0 * acc:.1f}%"


def _line(label: str, b: dict[str, Any]) -> str:
    n_en, n_zh, n_dual = b["n_en"], b["n_zh"], b["n_dual"]
    en = f"{b['ok_en']}/{n_en} ({_fmt_pct(b['acc_en'])})" if n_en else "-"
    zh = f"{b['ok_zh']}/{n_zh} ({_fmt_pct(b['acc_zh'])})" if n_zh else "-"
    dual = f"{b['ok_dual']}/{n_dual} ({_fmt_pct(b['acc_dual'])})" if n_dual else "-"
    return f"  {label:<12} EN={en:<16} ZH={zh:<16} dual={dual}"


def print_report(report: dict[str, Any]) -> None:
    total = {k: report[k] for k in ("ok_en", "n_en", "ok_zh", "n_zh", "ok_dual", "n_dual",
                                    "acc_en", "acc_zh", "acc_dual")}
    print(_line("total", total))
    for name, b in (report.get("by_source") or {}).items():
        print(_line(f"src={name}", b))
    for name, b in (report.get("by_tier") or {}).items():
        print(_line(f"tier={name}", b))
    named = report.get("results") or []
    n = len(named)
    print(f"  named coverage: en={report['named_en']}/{n}  zh={report['named_zh']}/{n}  "
          f"fails(dual)={len(report['fails'])}")


# ------------------------------------------------------------------------ main

def run_engine_benchmark(
    engine: str,
    data_path: str | Path,
    limit: int | None = None,
    workers: int | None = None,
    timeout: float = 3.0,
    metric: str = "strict",
) -> dict[str, Any]:
    """Score one engine on the gold data with a process pool."""
    rows = _load_rows(Path(data_path), limit)
    n_rows = len(rows)
    n_workers = max(1, workers if workers is not None else _default_workers())
    os.environ["NAMER_BENCH_ENGINE"] = engine
    os.environ["NAMER_BENCH_METRIC"] = metric
    t0 = time.perf_counter()
    print(f"engine={engine} metric={metric} start n={n_rows} workers={n_workers}", flush=True)
    if n_rows == 0:
        return _aggregate([])
    results = _run_pool(rows, n_workers, timeout, t0)
    report = _aggregate(results)
    report["engine"] = engine
    report["metric"] = metric
    report["elapsed_sec"] = round(time.perf_counter() - t0, 2)
    return report


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Score an arbitrary naming engine on the gold benchmark")
    p.add_argument("--engine", choices=["v2", "v3", "src"], default="v3",
                   help="v2/v3 = tools/namepredict-v{2,3}; src = active src namer (default: v3)")
    p.add_argument("--data", type=Path, default=Path("benchmarks/merged_benchmark.json"))
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--timeout", type=float, default=3.0, help="Per-row hard timeout (s)")
    p.add_argument("--metric", choices=["strict", "loose"], default="strict",
                   help="strict=逐字(仅去空格); loose=再剥括号/连字符/位次数字")
    p.add_argument("--json", action="store_true", help="Print full report as JSON")
    p.add_argument("--save-fails", type=Path, default=None,
                   help="Optional path to dump failed rows (gold + pred)")
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    report = run_engine_benchmark(args.engine, args.data, args.limit, args.workers, args.timeout,
                                  metric=args.metric)
    if args.json:
        slim = {k: v for k, v in report.items() if k != "results"}
        print(json.dumps(slim, ensure_ascii=False, indent=1))
    else:
        print_report(report)
        print(f"elapsed={report['elapsed_sec']}s", flush=True)
    if args.save_fails:
        out = [{
            "smiles": f.get("smiles"), "source": f.get("source"), "tier": f.get("tier"),
            "gold_en": f.get("english_name"), "gold_zh": f.get("chinese_name"),
            "pred_en": f.get("pred_en"), "pred_zh": f.get("pred_zh"),
            "en_ok": f.get("en_ok"), "zh_ok": f.get("zh_ok"),
        } for f in report["fails"]]
        args.save_fails.parent.mkdir(parents=True, exist_ok=True)
        args.save_fails.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"fails written: {args.save_fails} ({len(out)})", flush=True)


if __name__ == "__main__":
    main()
