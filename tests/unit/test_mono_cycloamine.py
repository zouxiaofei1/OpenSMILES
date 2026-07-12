# IUPAC: P-62.2.1 / P-22.1.1
# Layer: L2,L4,L5
"""Unsubstituted monocyclic primary monoamines: cyclohexanamine etc.

Parent = saturated monocarbocycle with one ring primary amine.
Unsubstituted: omit locant (cyclohexanamine not cyclohexan-1-amine).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive
    ("NC1CCCCC1", "cyclohexanamine", "环己胺"),
    ("NC1CCCC1", "cyclopentanamine", "环戊胺"),
    ("NC1CCC1", "cyclobutanamine", "环丁胺"),
    ("NC1CC1", "cyclopropanamine", "环丙胺"),
    ("NC1CCCCCC1", "cycloheptanamine", "环庚胺"),
    # negative
    ("CCN", "ethanamine", "乙胺"),
    ("CCCN", "propan-1-amine", "丙-1-胺"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_cycloamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_amine_name() -> None:
    r = SMILESNNamer().name("NC1CCCCC1")
    assert "nonan" not in normalize_en(r.en)
    assert normalize_en(r.en) != "hexanamine"
