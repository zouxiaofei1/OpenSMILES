# IUPAC: P-14.3.4 / P-62.2.1 / P-61.5
# Layer: L2,L4,L5
"""Multi-substituted benzene (cap 4) and benzene-1,n-diamine parents.

Simple ring subs (halo / alkyl / alkoxy / nitro / CF3) up to 4 on benzene.
Two primary amines on benzene → benzene-a,b-diamine (≤2 simple halo/methyl).
Gold chinese_name wins when it conflicts with pure system zh (anisole-style).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: 4-sub benzene
    (
        "ClC1=C(C(=CC(=C1)OC)[N+](=O)[O-])C",
        "1-chloro-5-methoxy-2-methyl-3-nitrobenzene",
        "1-氯-5-甲氧基-2-甲基-3-硝基苯",
    ),
    (
        "FC1=C(C=C(C(=C1)F)[N+](=O)[O-])[N+](=O)[O-]",
        "1,5-difluoro-2,4-dinitrobenzene",
        "1,5-二氟-2,4-二硝基苯",
    ),
    # negative: keep existing correct behaviour
    ("Clc1ccc(Cl)c(Cl)c1", "1,2,4-trichlorobenzene", "1,2,4-三氯苯"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzene_poly_diamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
