# IUPAC: P-29.3.1 / P-22.1.3
# Layer: L2,L3,L5
"""Simple mono-substituted alkylbenzenes: C1–C4 n-alkyl + isopropyl.

Parent = benzene (aromatic monocarbocycle). Monosubstituted: omit locant.
Methylbenzene retains toluene; isopropyl is branched alkyl (not propyl).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: mono C1–C4 n-alkyl + isopropylbenzene
    ("Cc1ccccc1", "toluene", "甲苯"),
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    ("CCCc1ccccc1", "propylbenzene", "丙基苯"),
    ("CCCCc1ccccc1", "butylbenzene", "丁基苯"),
    ("CC(C)c1ccccc1", "isopropylbenzene", "异丙基苯"),
    # negative: unsubstituted / monohalo / non-arene keep existing names
    ("c1ccccc1", "benzene", "苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)C", "2-methylpropane", "2-甲基丙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkylbenzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("CC(C)c1ccccc1", "2-methyloctane"),
        ("CCCc1ccccc1", "nonane"),
        ("CCCCc1ccccc1", "decane"),
    ],
)
def test_not_chain_alkane_name(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != normalize_en(bad)


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(C)Cc1ccccc1",   # isobutylbenzene
        "CCC(C)c1ccccc1",   # sec-butylbenzene
        "CC(C)(C)c1ccccc1", # tert-butylbenzene
    ],
)
def test_branched_butyl_not_bare_benzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert not (r.success and normalize_en(r.en) == "benzene")
