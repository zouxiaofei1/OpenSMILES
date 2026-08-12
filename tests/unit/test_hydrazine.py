# IUPAC: P-68.3.1.2
# Layer: L0,L1,L2,L5
"""Hydrazine scope (P-68.3.1.2) — negative guard only.

Positive hydrazine cases are not covered in this file. The retained cases
assert that acetamide and aniline are not named as hydrazines.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: amide / aniline / urea
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hydrazine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
