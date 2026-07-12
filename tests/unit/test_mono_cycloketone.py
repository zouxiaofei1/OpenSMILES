# IUPAC: P-64.2.1 / P-22.1.1
# Layer: L2,L4,L5
"""Unsubstituted monocyclic monoketones: cyclohexanone etc.

Parent = saturated monocarbocycle with one ring carbonyl.
Unsubstituted: omit locant (cyclohexanone not cyclohexan-1-one).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted cycloketones (omit locant)
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
    ("O=C1CCCC1", "cyclopentanone", "环戊酮"),
    ("O=C1CCC1", "cyclobutanone", "环丁酮"),
    ("O=C1CCCCCC1", "cycloheptanone", "环庚酮"),
    ("O=C1CC1", "cyclopropanone", "环丙酮"),
    # negative: chain ketones and other parents stay correct
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("CC(=O)CC", "butan-2-one", "丁-2-酮"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("CCO", "ethanol", "乙醇"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_cycloketone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_ketone_name() -> None:
    r = SMILESNNamer().name("O=C1CCCCC1")
    assert normalize_en(r.en) != "hexanone"
    assert "nonan" not in normalize_en(r.en)
