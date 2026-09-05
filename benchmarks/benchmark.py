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

from namepredict.constants import nospace, normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_MISSING_DATA_HINT = (
    "Generate with: python tools/merge_datasets.py "
    "--tiers smiles_tiers.json --chebi chebi20_test_1k.json "
    "--out data/merged_benchmark.json"
)

# Last-run snapshot for diff vs previous correct/status results.
_SNAPSHOT_PATH = Path(__file__).resolve().parent / ".last_run.json"


def _check_en(pred: str, gold: str) -> bool:
    # 空格不敏感：gold/预测在分词空格上可存在书写差异（如 acid 盐后缀前空格）
    return nospace(normalize_en(pred)) == nospace(normalize_en(gold))


def _check_zh(pred: str, gold: str) -> bool:
    # 空格不敏感：中文名内部无空格与空格并存（如 X盐酸盐 / X 盐酸盐）
    return nospace(normalize_zh(pred)) == nospace(normalize_zh(gold))


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


def _row_key(row: dict[str, Any]) -> str:
    rid = row.get("id")
    if rid is not None and str(rid) != "":
        return f"id:{rid}"
    return f"smiles:{row.get('smiles') or ''}"


def _result_entry(row: dict[str, Any], score: dict[str, Any], pred_en: str, pred_zh: str) -> dict:
    return {
        "key": _row_key(row),
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
    results: list[dict],
) -> None:
    score, pred_en, pred_zh = _score_one(namer, row)
    _tally(total, score)
    _tally(by_source[str(row.get("source") or "unknown")], score)
    _tally(by_tier[row.get("tier", 0)], score)
    entry = _result_entry(row, score, pred_en, pred_zh)
    results.append(entry)
    if _is_fail(score):
        fails.append(entry)


def _finalize(
    total: dict[str, int],
    by_source: dict[str, dict[str, int]],
    by_tier: dict[Any, dict[str, int]],
    fails: list[dict],
    results: list[dict] | None = None,
) -> dict[str, Any]:
    report = _bucket_report(total)
    report["fails"] = fails
    report["results"] = results or []
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
    results: list[dict] = []
    for row in rows:
        _process_row(namer, row, total, by_source, by_tier, fails, results)
    return _finalize(total, by_source, by_tier, fails, results)


def _fmt_pct(acc: float, ok: int, n: int) -> str:
    return f"{100.0 * acc:.1f}% ({ok}/{n})"


def _print_summary(report: dict[str, Any]) -> None:
    en_s = _fmt_pct(report["acc_en"], report["ok_en"], report["n_en"])
    zh_s = _fmt_pct(report["acc_zh"], report["ok_zh"], report["n_zh"])
    dual_s = _fmt_pct(report["acc_dual"], report["ok_dual"], report["n_dual"])
    print(f"en={en_s} zh={zh_s} dual={dual_s} fails={len(report['fails'])}")


def _snapshot_item(entry: dict[str, Any]) -> dict[str, Any]:
    return {
        "dual_ok": bool(entry.get("dual_ok")),
        "en_ok": entry.get("en_ok"),
        "zh_ok": entry.get("zh_ok"),
        "pred_en": entry.get("pred_en") or "",
        "pred_zh": entry.get("pred_zh") or "",
        "smiles": entry.get("smiles") or "",
        "english_name": entry.get("english_name") or "",
        "chinese_name": entry.get("chinese_name") or "",
    }


def _build_snapshot(report: dict[str, Any]) -> dict[str, Any]:
    items = {
        str(e.get("key") or _row_key(e)): _snapshot_item(e)
        for e in (report.get("results") or [])
    }
    return {
        "ok_dual": int(report.get("ok_dual") or 0),
        "n_dual": int(report.get("n_dual") or 0),
        "items": items,
    }


def _load_snapshot(path: Path = _SNAPSHOT_PATH) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) and isinstance(data.get("items"), dict) else None


def _save_snapshot(report: dict[str, Any], path: Path = _SNAPSHOT_PATH) -> None:
    snap = _build_snapshot(report)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(snap, f, ensure_ascii=False)


def _status_tuple(item: dict[str, Any]) -> tuple:
    return (item.get("dual_ok"), item.get("en_ok"), item.get("zh_ok"),
            item.get("pred_en") or "", item.get("pred_zh") or "")


def _diff_kind(prev: dict[str, Any], cur: dict[str, Any]) -> str | None:
    if _status_tuple(prev) == _status_tuple(cur):
        return None
    p_ok, c_ok = bool(prev.get("dual_ok")), bool(cur.get("dual_ok"))
    if p_ok and not c_ok:
        return "REGRESS"
    if not p_ok and c_ok:
        return "IMPROVE"
    return "CHANGE"


def _compare_to_previous(report: dict[str, Any], prev: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not prev:
        return []
    prev_items = prev.get("items") or {}
    diffs: list[dict[str, Any]] = []
    for entry in report.get("results") or []:
        key = str(entry.get("key") or _row_key(entry))
        old = prev_items.get(key)
        if not isinstance(old, dict):
            continue
        kind = _diff_kind(old, entry)
        if kind:
            diffs.append({"kind": kind, "key": key, "prev": old, "cur": entry})
    return diffs


def _fmt_ok(v: Any) -> str:
    if v is True:
        return "ok"
    if v is False:
        return "fail"
    return "-"


def _print_one_diff(d: dict[str, Any]) -> None:
    cur, prev = d["cur"], d["prev"]
    smiles = cur.get("smiles") or prev.get("smiles") or ""
    print(
        f"  [{d['kind']}] {d['key']}  smiles={smiles}  "
        f"dual {_fmt_ok(prev.get('dual_ok'))}->{_fmt_ok(cur.get('dual_ok'))}  "
        f"en {_fmt_ok(prev.get('en_ok'))}->{_fmt_ok(cur.get('en_ok'))}  "
        f"zh {_fmt_ok(prev.get('zh_ok'))}->{_fmt_ok(cur.get('zh_ok'))}"
    )
    print(f"    gold_en={cur.get('english_name') or prev.get('english_name') or ''!r}")
    print(f"    gold_zh={cur.get('chinese_name') or prev.get('chinese_name') or ''!r}")
    print(f"    prev_pred_en={prev.get('pred_en') or ''!r}  cur_pred_en={cur.get('pred_en') or ''!r}")
    print(f"    prev_pred_zh={prev.get('pred_zh') or ''!r}  cur_pred_zh={cur.get('pred_zh') or ''!r}")


_MAX_PRINT_DIFFS = 200


def _print_diffs(diffs: list[dict[str, Any]]) -> None:
    if not diffs:
        print("diff_vs_last: none (same as last run, or no previous snapshot)")
        return
    n_reg = sum(1 for d in diffs if d["kind"] == "REGRESS")
    n_imp = sum(1 for d in diffs if d["kind"] == "IMPROVE")
    n_chg = sum(1 for d in diffs if d["kind"] == "CHANGE")
    print(f"diff_vs_last: total={len(diffs)} REGRESS={n_reg} IMPROVE={n_imp} CHANGE={n_chg}")
    shown = diffs[:_MAX_PRINT_DIFFS]
    for d in shown:
        _print_one_diff(d)
    n_hidden = len(diffs) - len(shown)
    if n_hidden > 0:
        print(f"  ... {n_hidden} more diff(s) not shown (cap {_MAX_PRINT_DIFFS})")


def _report_for_json(report: dict[str, Any], diffs: list[dict[str, Any]]) -> dict[str, Any]:
    """Drop bulky per-row results; keep compact diff summary for --json consumers."""
    out = {k: v for k, v in report.items() if k != "results"}
    out["diff_vs_last"] = [
        {
            "kind": d["kind"],
            "key": d["key"],
            "smiles": (d["cur"].get("smiles") or d["prev"].get("smiles") or ""),
            "prev_dual_ok": d["prev"].get("dual_ok"),
            "cur_dual_ok": d["cur"].get("dual_ok"),
            "prev_pred_en": d["prev"].get("pred_en") or "",
            "cur_pred_en": d["cur"].get("pred_en") or "",
            "prev_pred_zh": d["prev"].get("pred_zh") or "",
            "cur_pred_zh": d["cur"].get("pred_zh") or "",
        }
        for d in diffs
    ]
    return out


def _handle_report(report: dict[str, Any], as_json: bool, snapshot_path: Path = _SNAPSHOT_PATH) -> None:
    prev = _load_snapshot(snapshot_path)
    diffs = _compare_to_previous(report, prev)
    report["diff_vs_last"] = diffs
    if as_json:
        print(json.dumps(_report_for_json(report, diffs), ensure_ascii=False))
    else:
        _print_summary(report)
        _print_diffs(diffs)
    _save_snapshot(report, snapshot_path)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Bilingual SMILES namer benchmark")
    p.add_argument("--data", type=Path, required=True, help="Path to merged_benchmark.json")
    p.add_argument("--limit", type=int, default=None, help="Optional row limit")
    p.add_argument("--json", action="store_true", help="Print full report as JSON")
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
    try:
        report = run_benchmark(SMILESNNamer(), args.data, limit=args.limit)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
    if args.no_snapshot:
        if args.json:
            print(json.dumps({k: v for k, v in report.items() if k != "results"}, ensure_ascii=False))
        else:
            _print_summary(report)
        return
    _handle_report(report, as_json=args.json, snapshot_path=args.snapshot)


if __name__ == "__main__":
    main()
