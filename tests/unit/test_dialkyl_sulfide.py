# IUPAC: P-63.2.1
# Layer: L1,L2,L5
"""Open-chain dialkyl sulfide scope (P-63.2.1) — negative guard only.

Positive sulfide cases are not covered in this file. The retained cases
assert ethanol / sulfoxide / thioester are not named as dialkyl sulfides.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: thiol / ether / alcohol must not break
    ("CCO", "ethanol", "乙醇"),
]

# sulfoxide / thioester must not be named as dialkyl sulfide
NEG_NOT_SULFIDE = [
    ("CS(C)=O", "dimethyl sulfide"),
    ("CC(=O)SC", "methyl sulfide"),
]


@pytest.mark.parametrize("smiles,forbidden_en", NEG_NOT_SULFIDE)
def test_not_dialkyl_sulfide(smiles: str, forbidden_en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert en != normalize_en(forbidden_en)
    assert "sulfide" not in en
    assert "硫醚" not in normalize_zh(r.zh)


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_dialkyl_sulfide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
