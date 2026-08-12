# IUPAC: P-66.5.1 / P-22.2.1
# Layer: L2,L3,L4,L5
"""Pyridinecarbonitrile scope (P-66.5.1) — negative guard only.

Positive pyridinecarbonitrile cases are not covered in this file. The retained
cases assert benzonitrile / pyridine / acetonitrile / propanenitrile stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: must not steal benzonitrile / pyridine / acid / chain nitrile
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CC#N", "acetonitrile", "乙腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinecarbonitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
