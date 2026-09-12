# IUPAC: P-64.1.1
# Layer: L2,L4,L5
"""Acetophenone scope (P-64.1.1) — negative guard only.

Positive acetophenone cases are not covered in this file. The retained cases
assert benzaldehyde / propan-2-one / benzene / benzoic acid / ethanol are not
named as acetophenone.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: benzaldehyde, acetone, benzene, benzoic acid, ethanol
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("c1ccccc1", "benzene", "苯"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_acetophenone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
