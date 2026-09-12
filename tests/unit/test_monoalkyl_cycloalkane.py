# IUPAC: P-14.3.4 / P-22.1.1
# Layer: L2,L3,L4,L5
"""Monosubstituted monocycloalkanes: methylcyclohexane etc.

Parent = saturated monocarbocycle; single linear alkyl side chain.
Monosubstituted: omit locant (methylcyclohexane not 1-methyl…).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: monoalkyl cycloalkanes (omit locant)
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("CC1CCCC1", "methylcyclopentane", "甲基环戊烷"),
    ("CC1CCC1", "methylcyclobutane", "甲基环丁烷"),
    ("CC1CC1", "methylcyclopropane", "甲基环丙烷"),
    ("CCC1CCCCC1", "ethylcyclohexane", "乙基环己烷"),
    # unsubstituted still ok
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("C1CC1", "cyclopropane", "环丙烷"),
    # negative
    ("CCC", "propane", "丙烷"),
    ("CCCCCC", "hexane", "己烷"),
    ("CC(C)C", "2-methylpropane", "2-甲基丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CCN", "ethanamine", "乙胺"),
    ("CC#N", "acetonitrile", "乙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_monoalkyl_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_dimethyl_not_methylcyclohexane() -> None:
    """Disubstituted ring out of simple monoalkyl scope naming as methylcyclohexane."""
    r = SMILESNNamer().name("CC1CCCCC1C")
    assert normalize_en(r.en) != "methylcyclohexane"
