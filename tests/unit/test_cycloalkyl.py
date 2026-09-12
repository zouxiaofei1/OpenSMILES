# IUPAC: P-29.6 / P-22.1.1 / P-14.3.4
# Layer: L2,L3,L5
"""Monocycloalkyl substituents (unsubstituted C3–C8) on cycloalkane parents.

P-29.6: cycloalkyl prefixes (cyclohexyl, cyclopentyl, …).
Multi-ring molecules pick one sat carbocycle as parent; the other is a side.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # P0: bi(cyclohexane)
    ("C1CCC(CC1)C1CCCCC1", "cyclohexylcyclohexane", "环己基环己烷"),
    # P1: 1-cyclohexylethyl
    (
        "C1(CCCCC1)C(C)C1CCCCC1",
        "(1-cyclohexylethyl)cyclohexane",
        "(1-环己基乙基)环己烷",
    ),
    (
        "CC(C1CCCCC1)C1CCCCC1",
        "(1-cyclohexylethyl)cyclohexane",
        "(1-环己基乙基)环己烷",
    ),
    # generality: cyclopentyl
    ("C1CCC(CC1)C1CCCC1", "cyclopentylcyclohexane", "环戊基环己烷"),
    # regressions
    ("C1CCC(CC1)C", "methylcyclohexane", "甲基环己烷"),
    ("CC(C)(C)C1CCCCC1", "tert-butylcyclohexane", "叔丁基环己烷"),
    (
        "C1CC(C(C)(CC)C)CCC1",
        "(2-methylbutan-2-yl)cyclohexane",
        "(2-甲基丁-2-基)环己烷",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkyl(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
