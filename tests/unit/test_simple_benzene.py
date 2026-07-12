# IUPAC: P-22.1.3
# Layer: L2,L4,L5
"""Retained name benzene: unsubstituted + monohalo + mono C1–C2 n-alkylbenzene.

Parent = aromatic monocarbocycle (benzene). Monosubstituted: omit locant.
Methylbenzene retains toluene; ethylbenzene is systematic.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted + monohalo + mono C1–C2 alkyl
    ("c1ccccc1", "benzene", "苯"),
    ("Fc1ccccc1", "fluorobenzene", "氟苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("Brc1ccccc1", "bromobenzene", "溴苯"),
    ("Ic1ccccc1", "iodobenzene", "碘苯"),
    ("Cc1ccccc1", "toluene", "甲苯"),
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    # negative: saturated cyclo / chain alcohol must stay correct
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_benzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_alkane_name() -> None:
    r = SMILESNNamer().name("c1ccccc1")
    assert normalize_en(r.en) != "hexane"


@pytest.mark.parametrize("smiles", ["C=Cc1ccccc1", "C#Cc1ccccc1"])
def test_unsaturated_sidechain_not_ethylbenzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != "ethylbenzene"
