# IUPAC: P-61.3.1 / P-14.3.4 / P-22.1.1
# Layer: L2,L3,L5
"""Monosubstituted monohalo monocycloalkanes: chlorocyclohexane etc.

Parent = saturated monocarbocycle; single halogen on ring carbon.
Monosubstituted: omit locant (chlorocyclohexane not 1-chlorocyclohexane).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: monohalo cycloalkanes (omit locant)
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
    ("BrC1CCCCC1", "bromocyclohexane", "溴环己烷"),
    ("FC1CCCCC1", "fluorocyclohexane", "氟环己烷"),
    ("IC1CCCCC1", "iodocyclohexane", "碘环己烷"),
    ("ClC1CCCC1", "chlorocyclopentane", "氯环戊烷"),
    ("ClC1CCC1", "chlorocyclobutane", "氯环丁烷"),
    ("ClC1CC1", "chlorocyclopropane", "氯环丙烷"),
    # unsubstituted / alkyl still ok
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    # negative: chain halo and other parents
    ("CCCCCCCl", "1-chlorohexane", "1-氯己烷"),
    ("ClCC", "chloroethane", "氯乙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_monohalo_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_halo_name() -> None:
    r = SMILESNNamer().name("ClC1CCCCC1")
    assert normalize_en(r.en) != "1-chlorohexane"
