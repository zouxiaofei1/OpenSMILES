# IUPAC: P-65.3
# Layer: L1,L2,L5
"""N-alkyl / N,N-dialkyl sulfonamide scope (P-65.3) — negative guard only.

Positive N-alkyl sulfonamide cases are not covered in this file. The retained
cases assert sulfoxide / amide / sulfone / sulfonate ester are not named as
N-alkyl sulfonamides.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer


# Must not interfere with existing
NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),       # sulfoxide
    ("CC(=O)N", "sulfonamide"),       # amide
    ("CS(=O)(=O)C", "sulfonamide"),   # sulfone
    ("CS(=O)(=O)OC", "sulfonamide"),  # sulfonate ester
]


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide_n_alkyl(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en
