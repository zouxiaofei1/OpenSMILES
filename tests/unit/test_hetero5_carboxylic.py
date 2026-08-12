# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
"""Hetero5 carboxylic acid scope (P-65.1.1) — negative guard only.

Positive hetero5-acid cases are not covered in this file. The retained case
asserts benzoic acid is not renamed as a hetero5 carboxylic acid.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: must not steal existing correct parents
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hetero5_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
