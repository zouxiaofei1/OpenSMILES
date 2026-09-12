# IUPAC: P-63.1.4 / P-62.5
# Layer: L2,L3,L4,L5
"""Aminophenols: phenol parent + amino prefix (OH senior to NH2).

Benzene with exactly one phenolic OH and one primary ring NH2.
OH is locant 1; amino gets 2/3/4. Optional extra halo/methyl/nitro
still limited by arene FG sub cap (amino counts as one extra).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted aminophenols (o/m/p)
    ("Nc1ccc(O)cc1", "4-aminophenol", "4-氨基苯酚"),
    ("Nc1ccccc1O", "2-aminophenol", None),
    ("Nc1cc(O)ccc1", "3-aminophenol", None),
    ("Oc1ccc(N)cc1", "4-aminophenol", "4-氨基苯酚"),  # isomorphic to first
    # negative near-miss: plain phenol / aniline / nitrophenol / chain / bare
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("Oc1ccc([N+]([O-])=O)cc1", "4-nitrophenol", "4-硝基苯酚"),
    ("CCO", "ethanol", "乙醇"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_aminophenol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
