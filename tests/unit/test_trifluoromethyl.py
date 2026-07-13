# IUPAC: P-29.3.2.2 / P-14.3.4
# Layer: L2,L3,L5
"""Aromatic trifluoromethyl (–CF₃) as a simple ring substituent prefix.

Ring-carbon-attached C with exactly three F and no further carbon is one
trifluoromethyl substituent (not three fluoro + alkyl). Allowed on retained
arene parents (benzene, phenol, benzoic, …) within the existing ≤3 sub cap.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: mono CF3 benzene
    ("FC(F)(F)c1ccccc1", "trifluoromethylbenzene", "三氟甲基苯"),
    # positive: phenol + halo + CF3
    ("FC=1C=CC(=C(C1)O)C(F)(F)F", "5-fluoro-2-(trifluoromethyl)phenol", "5-氟-2-三氟甲基苯酚"),
    # positive: hydroxy-benzoic + CF3
    (
        "OC=1C=C(C(=O)O)C=C(C1)C(F)(F)F",
        "3-hydroxy-5-(trifluoromethyl)benzoic acid",
        "3-羟基-5-三氟甲基苯甲酸",
    ),
    # positive: tri-sub benzene halo + CF3
    (
        "FC1=C(C(=CC=C1)C(F)(F)F)I",
        "1-fluoro-2-iodo-3-(trifluoromethyl)benzene",
        "1-氟-2-碘-3-三氟甲基苯",
    ),
    # negative: keep existing correct names / chain CF3 not arene rule
    ("Cc1ccccc1", "toluene", "甲苯"),
    ("Fc1ccccc1", "fluorobenzene", "氟苯"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CC(F)(F)F", "1,1,1-trifluoroethane", "1,1,1-三氟乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_trifluoromethyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
