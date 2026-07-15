# IUPAC: P-65.6 / carbonate
# Layer: L1,L2,L5
"""Organic carbonate RO–C(=O)–OR′ functional-class names (P-65.6).

Scope (first cut):
- symmetric dialkyl R=R′ = n-alkyl C1–C4: dimethyl / diethyl carbonate
- alkyl–aryl: R = n-alkyl C1–C4, Ar = unfused Ph with ≤2 Me/halo

Out of scope: inorganic carbonate salts, cyclic carbonates, dipyridyl carbonate.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # symmetric dialkyl
    ("COC(=O)OC", "dimethyl carbonate", "碳酸二甲酯"),
    ("CCOC(=O)OCC", "diethyl carbonate", "碳酸二乙酯"),
    # alkyl–aryl (bench gold ZH: aryl + alkyl + 碳酸酯)
    (
        "C(OC)(OC1=C(C=C(C=C1)C)C)=O",
        "methyl 2,4-dimethylphenyl carbonate",
        "2,4-二甲苯基甲基碳酸酯",
    ),
    # near-miss negatives: keep existing correct names
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("CCOC(=O)N", "ethyl carbamate", "氨基甲酸乙酯"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
]


# Must not be named as carbonate
NEG_NOT_CARBONATE = [
    ("CCOC(=O)C", "carbonate"),
    ("CCOC(=O)N", "carbonate"),
    ("CC(=O)OC", "carbonate"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_organic_carbonate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_CARBONATE)
def test_not_carbonate(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "碳酸" not in zh or "氨基甲酸" in zh
