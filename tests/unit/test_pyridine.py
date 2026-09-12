# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained pyridine scope (P-22.2.1) — unsubstituted + negative guard.

Only the unsubstituted pyridine positive case is retained here; substituted /
carboxylic cases are not covered. The other cases assert benzene / ethanol /
cyclohexane stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted + monomethyl + monohalo
    ("c1ccncc1", "pyridine", "吡啶"),
    # negative: carbocycles / simple chain must stay correct
    ("c1ccccc1", "benzene", "苯"),
    ("CCO", "ethanol", "乙醇"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_pyridine_not_pentane() -> None:
    r = SMILESNNamer().name("c1ccncc1")
    assert r.success
    assert normalize_en(r.en) != "pentane"
