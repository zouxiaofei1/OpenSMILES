# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
"""Sat-monohetero carboxylic acid scope (P-65.1.1) — negative guard only.

Positive sat-hetero acid cases are not covered in this file. The retained cases
assert benzoic acid and hexanoic acid are not renamed as sat-hetero acids.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: must not steal existing correct parents
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("O=C(O)CCCCC", "hexanoic acid", "己酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
