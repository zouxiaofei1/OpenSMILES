# IUPAC: P-56.4 / P-29 (H2C= -> methylidene, preferred prefix; 不用旧名 methylene)
# Layer: L3,L5
"""Exocyclic double bond leaves: parent=CH2 must keep the double bond (methylidene).

A substituent leaf whose root atom bonds to the parent via a DOUBLE bond was
previously collapsed to the saturated alkyl (methyl) because the anchor dummy
bond was hardcoded single (methylcyclohexane for C=C1CCCCC1, formula loses H2).
The exocyclic =CH2 group is IUPAC 'methylidene' (P-56.4); endocyclic/chain
alkenes and saturated twins must stay unchanged.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

POSITIVE = [
    # exocyclic =CH2 on a ring -> methylidene retained
    ("C=C1CCCCC1", "methylidenecyclohexane", "亚甲基环己烷"),
    # exocyclic =CH2 on a chain near FG + stereo
    ("C=C(C[C@H](N)C(=O)O)C(=O)O",
     "(2S)-2-amino-4-methylidenepentanedioic acid",
     "(2S)-2-氨基-4-亚甲基戊二酸"),
    # nested: =CH2 on a cyclopropyl substituent (recursive path)
    ("C=C1CC1CC(=O)C(=O)O",
     "3-(2-methylidenecyclopropyl)-2-oxopropanoic acid",
     "3-(2-亚甲基环丙基)-2-氧代丙酸"),
]

NEGATIVE = [
    # saturated twin must NOT take the ylidene route
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    # endocyclic/chain alkenes stay alkene names
    ("C1=CCCCC1", "cyclohexene", "环己烯"),
    ("C1=CCCC1", "cyclopentene", "环戊烯"),
    ("C=C(C)C", "2-methylprop-1-ene", "2-甲基丙-1-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", POSITIVE)
def test_exocyclic_methylidene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", NEGATIVE)
def test_neighbors_not_ylidene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_exocyclic_not_saturated() -> None:
    """Guard: exocyclic =CH2 must not be reported as the saturated methyl twin."""
    r = SMILESNNamer().name("C=C1CCCCC1")
    assert normalize_en(r.en) != normalize_en("methylcyclohexane")
    assert "methylidene" in normalize_en(r.en)
