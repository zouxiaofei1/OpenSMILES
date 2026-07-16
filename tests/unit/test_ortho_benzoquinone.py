# IUPAC: P-64.2
# Layer: L2,L4,L5
"""1,2-Benzoquinone PIN as cyclohexa-3,5-diene-1,2-dione (P-64.2 Round C).

Single C6 carbocycle + exactly two ortho ring ketones + two endocyclic C=C.
PIN is systematic (not retained o-benzoquinone). First-cut: halo / OH / Me
or n-alkyl C1–C12 / alkoxy C1–C2. No naphthoquinone this round.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # unsubstituted parent (was 1,4-BQ negative methanone)
    (
        "O=C1C=CC=CC1=O",
        "cyclohexa-3,5-diene-1,2-dione",
        "环己-3,5-二烯-1,2-二酮",
    ),
    # gold: 4-hydroxy-5-methyl
    (
        "CC1=CC(=O)C(=O)C=C1O",
        "4-hydroxy-5-methylcyclohexa-3,5-diene-1,2-dione",
        "4-羟基-5-甲基环己-3,5-二烯-1,2-二酮",
    ),
    # simple methyl
    (
        "CC1=CC=CC(=O)C1=O",
        "3-methylcyclohexa-3,5-diene-1,2-dione",
        "3-甲基环己-3,5-二烯-1,2-二酮",
    ),
    # 4,5-dichloro
    (
        "ClC1=CC(=O)C(=O)C=C1Cl",
        "4,5-dichlorocyclohexa-3,5-diene-1,2-dione",
        "4,5-二氯环己-3,5-二烯-1,2-二酮",
    ),
    # negatives: 1,4-BQ must stay; anthraquinone / acetophenone unchanged
    (
        "O=C1C=CC(=O)C=C1",
        "cyclohexa-2,5-diene-1,4-dione",
        "环己-2,5-二烯-1,4-二酮",
    ),
    ("CC(=O)c1ccccc1", "acetophenone", "苯乙酮"),
    ("O=C1c2ccccc2C(=O)c2ccccc12", "9,10-anthraquinone", "蒽醌"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ortho_benzoquinone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
