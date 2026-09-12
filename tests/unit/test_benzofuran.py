# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained benzofuran scope (P-22.2.1) — negative guard only.

Positive benzofuran cases are not covered in this file. The retained cases
assert indole / naphthalene / benzene must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: near neighbors must not regress
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzofuran_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
