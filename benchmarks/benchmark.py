"""Bilingual field-wise benchmark scorer (read-only gold data)."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

_SRC = Path(__file__).resolve().parents[1] / "src"
_src_str = str(_SRC)
if _src_str not in sys.path:
    sys.path.insert(0, _src_str)

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_MISSING_DATA_HINT = (
    "Generate with: python tools/merge_datasets.py "
    "--tiers smiles_tiers.json --chebi chebi20_test_1k.json "
    "--out data/merged_benchmark.json"
)


def _check_en(pred: str, gold: str) -> bool:
    return normalize_en(pred) == normalize_en(gold)


def _check_zh(pred: str, gold: str) -> bool:
    return normalize_zh(pred) == normalize_zh(gold)


def score_record(pred_en: str, pred_zh: str, row: dict[str, Any]) -> dict[str, Any]:
    """Score one prediction against a gold row (field-wise)."""
    en_ok = _check_en(pred_en, row.get("english_name") or "") if row.get("eval_en") else None
    zh_ok = _check_zh(pred_zh, row.get("chinese_name") or "") if row.get("eval_zh") else None
    checks = [x for x in (en_ok, zh_ok) if x is not None]
    dual_ok = all(checks) if checks else False
    return {"en_ok": en_ok, "zh_ok": zh_ok, "dual_ok": dual_ok}


def _load_rows(data_path: Path, limit: int | None) -> list[dict[str, Any]]:
    if not data_path.is_file():
        raise FileNotFoundError(f"Benchmark data not found: {data_path}\n{_MISSING_DATA_HINT}")
    with data_path.open(encoding="utf-8") as f:
        rows = json.load(f)
    if not isinstance(rows, list):
        raise ValueError(f"Expected JSON list in {data_path}")
    return rows[:limit] if limit is not None else rows


def _empty_bucket() -> dict[str, int]:
    return {"ok_en": 0, "n_en": 0, "ok_zh": 0, "n_zh": 0, "ok_dual": 0, "n_dual": 0}


def _acc(ok: int, n: int) -> float:
    return 0.0 if n == 0 else ok / n


def _bucket_report(b: dict[str, int]) -> dict[str, float | int]:
    return {
        "acc_en": _acc(b["ok_en"], b["n_en"]),
        "acc_zh": _acc(b["ok_zh"], b["n_zh"]),
        "acc_dual": _acc(b["ok_dual"], b["n_dual"]),
        "n_en": b["n_en"], "n_zh": b["n_zh"], "n_dual": b["n_dual"],
        "ok_en": b["ok_en"], "ok_zh": b["ok_zh"], "ok_dual": b["ok_dual"],
    }


def _tally(bucket: dict[str, int], score: dict[str, Any]) -> None:
    if score["en_ok"] is not None:
        bucket["n_en"] += 1
        bucket["ok_en"] += int(score["en_ok"])
    if score["zh_ok"] is not None:
        bucket["n_zh"] += 1
        bucket["ok_zh"] += int(score["zh_ok"])
    if score["en_ok"] is not None or score["zh_ok"] is not None:
        bucket["n_dual"] += 1
        bucket["ok_dual"] += int(score["dual_ok"])


def _fail_entry(row: dict[str, Any], score: dict[str, Any], pred_en: str, pred_zh: str) -> dict:
    return {
        "id": row.get("id"), "smiles": row.get("smiles"),
        "source": row.get("source"), "tier": row.get("tier"),
        "features": list(row.get("features") or []),
        "english_name": row.get("english_name"), "chinese_name": row.get("chinese_name"),
        "pred_en": pred_en, "pred_zh": pred_zh,
        "en_ok": score["en_ok"], "zh_ok": score["zh_ok"], "dual_ok": score["dual_ok"],
    }


def _score_one(namer: SMILESNNamer, row: dict[str, Any]) -> tuple[dict, str, str]:
    result = namer.name(str(row.get("smiles") or ""))
    pred_en, pred_zh = result.en or "", result.zh or ""
    return score_record(pred_en, pred_zh, row), pred_en, pred_zh


def _is_fail(score: dict[str, Any]) -> bool:
    if score["en_ok"] is None and score["zh_ok"] is None:
        return False
    return not score["dual_ok"]


def _process_row(
    namer: SMILESNNamer,
    row: dict[str, Any],
    total: dict[str, int],
    by_source: dict[str, dict[str, int]],
    by_tier: dict[Any, dict[str, int]],
    fails: list[dict],
) -> None:
    score, pred_en, pred_zh = _score_one(namer, row)
    _tally(total, score)
    _tally(by_source[str(row.get("source") or "unknown")], score)
    _tally(by_tier[row.get("tier", 0)], score)
    if _is_fail(score):
        fails.append(_fail_entry(row, score, pred_en, pred_zh))


def _finalize(
    total: dict[str, int],
    by_source: dict[str, dict[str, int]],
    by_tier: dict[Any, dict[str, int]],
    fails: list[dict],
) -> dict[str, Any]:
    report = _bucket_report(total)
    report["fails"] = fails
    report["by_source"] = {k: _bucket_report(v) for k, v in sorted(by_source.items())}
    report["by_tier"] = {
        str(k): _bucket_report(v) for k, v in sorted(by_tier.items(), key=lambda x: str(x[0]))
    }
    return report


def run_benchmark(namer: SMILESNNamer, data_path: str | Path, limit: int | None = None) -> dict[str, Any]:
    """Run namer on gold JSON; return report with acc_* / n_* / fails / by_*."""
    rows = _load_rows(Path(data_path), limit)
    total = _empty_bucket()
    by_source: dict[str, dict[str, int]] = defaultdict(_empty_bucket)
    by_tier: dict[Any, dict[str, int]] = defaultdict(_empty_bucket)
    fails: list[dict] = []
    for row in rows:
        _process_row(namer, row, total, by_source, by_tier, fails)
    return _finalize(total, by_source, by_tier, fails)


def _fmt_pct(acc: float, ok: int, n: int) -> str:
    return f"{100.0 * acc:.1f}% ({ok}/{n})"


def _print_summary(report: dict[str, Any]) -> None:
    en_s = _fmt_pct(report["acc_en"], report["ok_en"], report["n_en"])
    zh_s = _fmt_pct(report["acc_zh"], report["ok_zh"], report["n_zh"])
    dual_s = _fmt_pct(report["acc_dual"], report["ok_dual"], report["n_dual"])
    print(f"en={en_s} zh={zh_s} dual={dual_s} fails={len(report['fails'])}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Bilingual SMILES namer benchmark")
    p.add_argument("--data", type=Path, required=True, help="Path to merged_benchmark.json")
    p.add_argument("--limit", type=int, default=None, help="Optional row limit")
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    try:
        report = run_benchmark(SMILESNNamer(), args.data, limit=args.limit)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
    _print_summary(report)


if __name__ == "__main__":
    main()
