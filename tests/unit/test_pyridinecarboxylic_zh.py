# IUPAC: P-65.1.1
# Layer: L5
"""Pyridinecarboxylic Chinese suffix scope (P-65.1.1) — negative guard only.

Positive 吡啶-n-甲酸 cases are not covered in this file. The retained case
asserts benzoic acid keeps 苯甲酸 (no regression to a pyridine suffix).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: carbocyclic carboxylic acids keep correct 甲酸 (no regression)
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinecarboxylic_zh(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
