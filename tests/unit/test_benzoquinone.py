# IUPAC: P-64.2.2.2.3
# Layer: L2,L4,L5
"""1,4-Benzoquinone PIN as cyclohexa-2,5-diene-1,4-dione (P-64.2).

Single C6 carbocycle + exactly two para ring ketones + two endocyclic
double bonds. PIN is systematic (not retained 1,4-benzoquinone).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # unsubstituted parent
    (
        "O=C1C=CC(=O)C=C1",
        "cyclohexa-2,5-diene-1,4-dione",
        "环己-2,5-二烯-1,4-二酮",
    ),
    # tetrachloro (chloranil systematic)
    (
        "O=C1C(Cl)=C(Cl)C(=O)C(Cl)=C1Cl",
        "2,3,5,6-tetrachlorocyclohexa-2,5-diene-1,4-dione",
        "2,3,5,6-四氯环己-2,5-二烯-1,4-二酮",
    ),
    # tetrahydroxy
    (
        "O=C1C(O)=C(O)C(=O)C(O)=C1O",
        "2,3,5,6-tetrahydroxycyclohexa-2,5-diene-1,4-dione",
        "2,3,5,6-四羟基环己-2,5-二烯-1,4-二酮",
    ),
    # hydroxy + methoxy + methyl (gold EN)
    (
        "COC1=C(O)C(=O)C(C)=CC1=O",
        "3-hydroxy-2-methoxy-5-methylcyclohexa-2,5-diene-1,4-dione",
        "3-羟基-2-甲氧基-5-甲基环己-2,5-二烯-1,4-二酮",
    ),
    # alt Kekulé / ring order for parent
    (
        "C1=CC(=O)C=CC1=O",
        "cyclohexa-2,5-diene-1,4-dione",
        "环己-2,5-二烯-1,4-二酮",
    ),
    # negatives: must not break existing correct names
    ("CC(=O)c1ccccc1", "acetophenone", "苯乙酮"),
    ("O=C1c2ccccc2C(=O)c2ccccc12", "9,10-anthraquinone", "蒽醌"),
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
    ("c1ccccc1", "benzene", "苯"),
    # open-chain dione stays alkanedione
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzoquinone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
