# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Dialkyl sulfone functional parents (P-65.3.1.2).

Scope:
- Symmetric dialkyl sulfone R-SO2-R (dimethyl sulfone, etc.)
- Unsymmetric dialkyl sulfone R-SO2-R'
- Sulfone as principal group (rank ~ ketone, 6)
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # dimethyl sulfone (symmetric)
    ("CS(=O)(=O)C", "dimethyl sulfone", "二甲基砜"),
    # ethyl methyl sulfone (unsymmetric)
    ("CCS(=O)(=O)C", "ethyl methyl sulfone", "乙基甲基砜"),
    # diethyl sulfone
    ("CCS(=O)(=O)CC", "diethyl sulfone", "二乙基砜"),
    # dipropyl sulfone
    ("CCCS(=O)(=O)CCC", "dipropyl sulfone", "二丙基砜"),
]

# Must not mis-identify
NEG_NOT_SULFONE = [
    ("CS(C)=O", "sulfone"),           # sulfoxide
    ("CS(=O)(=O)N", "sulfone"),       # sulfonamide
    ("CS(=O)(=O)O", "sulfone"),       # sulfonic acid
    ("CS(=O)(=O)Cl", "sulfone"),      # sulfonyl chloride
    ("CC(=O)C", "sulfone"),           # ketone (acetone)
    ("CS(=O)(=O)OC", "sulfone"),      # sulfonate ester
    ("CSC", "sulfone"),               # sulfide
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sulfone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"Failed for {smiles}: {r.en}"
    assert normalize_en(r.en) == normalize_en(en), f"{smiles}: got {r.en!r}, expected {en!r}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONE)
def test_not_sulfone(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en, f"{smiles} should not contain '{forbidden}' but got {r.en!r}"
