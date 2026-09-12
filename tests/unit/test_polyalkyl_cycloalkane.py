# IUPAC: P-14.3.4 / P-22.1.1
# Layer: L2,L4,L5
"""Polyalkyl monocyclic cycloalkanes: ≥2 linear C1–C4 alkyls on one ring.

Parent = saturated monocarbocycle; multiple linear alkyl side chains.
Lowest set of locants (P-14.3.4); multiplicative prefixes (P-29.2).
Monosubstituted still omits locant; unsubstituted and open-chain intact.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: polyalkyl monocycloalkanes (keep locants + mult prefixes)
    ("CC1CCCCC1C", "1,2-dimethylcyclohexane", "1,2-二甲基环己烷"),
    ("CC1CCC(C)CC1", "1,4-dimethylcyclohexane", "1,4-二甲基环己烷"),
    ("CC1CCCC(C)C1", "1,3-dimethylcyclohexane", "1,3-二甲基环己烷"),
    ("CC1(C)CCCCC1", "1,1-dimethylcyclohexane", "1,1-二甲基环己烷"),
    ("CC1CC(C)CC(C)C1", "1,3,5-trimethylcyclohexane", "1,3,5-三甲基环己烷"),
    ("CCC1CCCCC1C", "1-ethyl-2-methylcyclohexane", "1-乙基-2-甲基环己烷"),
    ("CC1CCCC1C", "1,2-dimethylcyclopentane", "1,2-二甲基环戊烷"),
    ("CCC1CCCCC1CC", "1,2-diethylcyclohexane", "1,2-二乙基环己烷"),
    # negative: mono / unsubstituted / open-chain / mono FG must not break
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("CCC1CCCCC1", "ethylcyclohexane", "乙基环己烷"),
    ("CCCC", "butane", "丁烷"),
    ("CCO", "ethanol", "乙醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyalkyl_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
