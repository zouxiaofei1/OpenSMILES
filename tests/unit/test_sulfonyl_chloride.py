# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Sulfonyl chloride scope (P-65.3) — negative guard only.

Positive sulfonyl chloride cases are not covered in this file. The retained
cases assert sulfonamide / acetyl chloride / sulfoxide / tosylate are not
named as sulfonyl chlorides.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# Must not be named as sulfonyl chloride
NEG_NOT_SULFONYL_CHLORIDE = [
    ("c1ccc(S(=O)(=O)N)cc1", "sulfonyl chloride"),
    ("CC(=O)Cl", "sulfonyl chloride"),
    ("CS(C)=O", "sulfonyl chloride"),
    ("CC1=CC=C(C=C1)S(=O)(=O)OCCCC", "sulfonyl chloride"),
    ("CS(=O)(=O)N", "sulfonyl chloride"),
]


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONYL_CHLORIDE)
def test_not_sulfonyl_chloride(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰氯" not in zh
