# IUPAC: P-65.3
# Layer: L1,L2,L5
"""N-alkyl and N,N-dialkyl sulfonamide functional parents (P-65.3).

Scope:
- N-alkyl alkanesulfonamide (secondary sulfonamide)
- N,N-dialkyl alkanesulfonamide (tertiary sulfonamide)
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # N-ethylmethanesulfonamide (secondary sulfonamide)
    ("CS(=O)(=O)NCC", "N-ethylmethanesulfonamide", "N-乙基甲磺酰胺"),
    # N,N-dimethylmethanesulfonamide (tertiary sulfonamide)
    ("CS(=O)(=O)N(C)C", "N,N-dimethylmethanesulfonamide", "N,N-二甲基甲磺酰胺"),
    # N-propylethanesulfonamide
    ("CCS(=O)(=O)NCCC", "N-propylethanesulfonamide", "N-丙基乙磺酰胺"),
    # N-butylpropanesulfonamide
    ("CCCS(=O)(=O)NCCCC", "N-butylpropanesulfonamide", "N-丁基丙磺酰胺"),
    # unsubstituted still works
    ("CS(=O)(=O)N", "methanesulfonamide", "甲磺酰胺"),
]

# Must not interfere with existing
NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),       # sulfoxide
    ("CC(=O)N", "sulfonamide"),       # amide
    ("CS(=O)(=O)C", "sulfonamide"),   # sulfone
    ("CS(=O)(=O)OC", "sulfonamide"),  # sulfonate ester
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_alkyl_sulfonamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide_n_alkyl(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en
