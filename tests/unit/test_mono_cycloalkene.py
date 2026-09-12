# IUPAC: P-31.1 / P-22.1.1
# Layer: L2,L4,L5
"""Unsubstituted monocyclic monoalkenes: cyclohexene etc.

Parent = monocarbocycle with one endocyclic C=C; no substituents.
Unsubstituted: omit locant (cyclohexene not cyclohex-1-ene).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive
    ("C1=CCCCC1", "cyclohexene", "环己烯"),
    ("C1=CCCC1", "cyclopentene", "环戊烯"),
    ("C1=CCC1", "cyclobutene", "环丁烯"),
    ("C1=CCCCCC1", "cycloheptene", "环庚烯"),
    # negative: chain alkenes and saturated rings
    ("C=CCC", "but-1-ene", "丁-1-烯"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("NC1CCCCC1", "cyclohexanamine", "环己胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_cycloalkene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_ene_name() -> None:
    r = SMILESNNamer().name("C1=CCCCC1")
    assert "dec" not in normalize_en(r.en)
    assert normalize_en(r.en) != "hexene"
