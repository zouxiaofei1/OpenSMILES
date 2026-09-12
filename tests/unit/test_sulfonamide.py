# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Sulfonamide scope (P-65.3) — negative guard only.

Positive sulfonamide cases are not covered in this file. The retained cases
assert acetamide / aniline / sulfoxide / sulfone are not named as sulfonamides.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative near-miss: must not become sulfonamide / keep correct names
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
]


# Must not be named as sulfonamide
NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),
    ("CC(=O)N", "sulfonamide"),
    ("c1ccccc1N", "sulfonamide"),
    ("CS(=O)(=O)C", "sulfonamide"),  # sulfone
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_sulfonamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰胺" not in zh
