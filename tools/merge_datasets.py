"""Merge smiles_tiers + chebi into a unified bilingual benchmark dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from rdkit import Chem


def canonical_smiles(smiles: str) -> str:
    """Return RDKit canonical SMILES, or the original string on failure."""
    if not smiles or not str(smiles).strip():
        return str(smiles or "")
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return str(smiles).strip()
    return Chem.MolToSmiles(mol)


def _to_record(raw: dict[str, Any], source: str) -> dict[str, Any]:
    english = str(raw.get("english_name") or "")
    chinese = str(raw.get("chinese_name") or "")
    prefix = "tiers" if source == "smiles_tiers" else "chebi"
    return {
        "id": f"{prefix}-{raw.get('id', '')}", "smiles": str(raw.get("smiles") or ""),
        "english_name": english, "chinese_name": chinese, "source": source,
        "tier": raw.get("tier", 1), "features": list(raw.get("features") or []),
        "eval_en": bool(english.strip()), "eval_zh": bool(chinese.strip()),
    }


def _should_replace(existing: dict[str, Any], new: dict[str, Any]) -> bool:
    e_zh = bool((existing.get("chinese_name") or "").strip())
    n_zh = bool((new.get("chinese_name") or "").strip())
    if n_zh and not e_zh:
        return True
    if e_zh and not n_zh:
        return False
    if new.get("source") == "smiles_tiers" and existing.get("source") != "smiles_tiers":
        return True
    return False


def _ingest(by_key: dict[str, dict[str, Any]], rows: list[dict[str, Any]], source: str) -> None:
    for raw in rows:
        rec = _to_record(raw, source)
        key = canonical_smiles(rec["smiles"])
        if key not in by_key or _should_replace(by_key[key], rec):
            by_key[key] = rec


def merge_records(
    tiers: list[dict[str, Any]],
    chebi: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Merge tiers + chebi lists keyed by canonical SMILES."""
    by_key: dict[str, dict[str, Any]] = {}
    _ingest(by_key, tiers, "smiles_tiers")
    _ingest(by_key, chebi, "chebi")
    return list(by_key.values())


def _load_json(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list):
        raise ValueError(f"Expected JSON list in {path}")
    return data


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Merge tiers + chebi benchmark datasets")
    parser.add_argument("--tiers", type=Path, required=True, help="Path to smiles_tiers.json")
    parser.add_argument("--chebi", type=Path, required=True, help="Path to chebi JSON")
    parser.add_argument("--out", type=Path, required=True, help="Output merged JSON path")
    return parser


def _write_and_report(merged: list[dict[str, Any]], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)
    n_zh = sum(1 for x in merged if x["eval_zh"])
    n_en = sum(1 for x in merged if x["eval_en"])
    print(f"Wrote {len(merged)} records to {out} (eval_zh={n_zh}, eval_en={n_en})")


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    tiers = _load_json(args.tiers)
    chebi = _load_json(args.chebi)
    merged = merge_records(tiers, chebi)
    _write_and_report(merged, args.out)


if __name__ == "__main__":
    main()
