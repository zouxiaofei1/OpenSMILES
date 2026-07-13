# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained name pyridine (azabenzene): unsubstituted, monomethyl, monohalo,
and pyridinecarboxylic acids with N = locant 1.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted + monomethyl + monohalo
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Cc1ccccn1", "2-methylpyridine", "2-甲基吡啶"),
    ("Clc1ccccn1", "2-chloropyridine", "2-氯吡啶"),
    # positive: pyridinecarboxylic acids (N=1)
    ("OC(=O)c1ccccn1", "pyridine-2-carboxylic acid", "吡啶-2-羧酸"),
    ("OC(=O)c1cccnc1", "pyridine-3-carboxylic acid", "吡啶-3-羧酸"),
    # negative: carbocycles / simple chain must stay correct
    ("c1ccccc1", "benzene", "苯"),
    ("Cc1ccccc1", "toluene", "甲苯"),
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


def test_methyl_locant_is_2_not_6() -> None:
    """N fixed as 1; monomethyl must be 2- not 6-."""
    r = SMILESNNamer().name("Cc1ccccn1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylpyridine"
    assert "6-methyl" not in en


def test_pyridine_not_pentane() -> None:
    r = SMILESNNamer().name("c1ccncc1")
    assert r.success
    assert normalize_en(r.en) != "pentane"
