# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Dialkyl sulfone scope (P-65.3.1.2) — negative guard only.

Positive sulfone cases are not covered in this file. The retained cases
assert sulfoxide / sulfonamide / sulfonic acid / sulfonyl chloride / ketone /
sulfonate ester / sulfide are not named as sulfones.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer


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


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONE)
def test_not_sulfone(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en, f"{smiles} should not contain '{forbidden}' but got {r.en!r}"
