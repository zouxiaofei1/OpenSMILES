from pathlib import Path

from tools.merge_datasets import merge_records


def test_merge_prefers_zh_and_sets_flags(tmp_path: Path):
    tiers = [
        {
            "id": 1,
            "smiles": "CCO",
            "english_name": "ethanol",
            "chinese_name": "乙醇",
            "tier": 1,
            "features": ["alcohol"],
        }
    ]
    chebi = [
        {
            "id": 9,
            "smiles": "CCO",
            "english_name": "ethanol",
            "chinese_name": "",
            "tier": 1,
            "features": [],
        }
    ]
    out = merge_records(tiers, chebi)
    assert len(out) == 1
    assert out[0]["eval_zh"] is True
    assert out[0]["chinese_name"] == "乙醇"
    assert out[0]["source"] == "smiles_tiers"


def test_merge_keeps_unique_and_chebi_only():
    tiers = [
        {
            "id": 1,
            "smiles": "C",
            "english_name": "methane",
            "chinese_name": "甲烷",
            "tier": 1,
            "features": [],
        }
    ]
    chebi = [
        {
            "id": 2,
            "smiles": "CC",
            "english_name": "ethane",
            "chinese_name": "",
            "tier": 1,
            "features": [],
        }
    ]
    out = merge_records(tiers, chebi)
    assert len(out) == 2
    by_smiles = {r["smiles"]: r for r in out}
    assert by_smiles["C"]["source"] == "smiles_tiers"
    assert by_smiles["C"]["id"] == "tiers-1"
    assert by_smiles["CC"]["source"] == "chebi"
    assert by_smiles["CC"]["id"] == "chebi-2"
    assert by_smiles["CC"]["eval_zh"] is False
    assert by_smiles["CC"]["eval_en"] is True


def test_merge_prefers_zh_from_chebi_over_empty_tiers():
    tiers = [
        {
            "id": 1,
            "smiles": "CCO",
            "english_name": "ethanol",
            "chinese_name": "",
            "tier": 1,
            "features": [],
        }
    ]
    chebi = [
        {
            "id": 9,
            "smiles": "CCO",
            "english_name": "ethanol",
            "chinese_name": "乙醇",
            "tier": 1,
            "features": [],
        }
    ]
    out = merge_records(tiers, chebi)
    assert len(out) == 1
    assert out[0]["chinese_name"] == "乙醇"
    assert out[0]["eval_zh"] is True
    assert out[0]["source"] == "chebi"
