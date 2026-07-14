# IUPAC: P-25
# Layer: L2,L5
"""Unsubstituted anthracene retained parent via ring_systems linear 666."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("c1ccc2cc3ccccc3cc2c1", "anthracene", "蒽"),
    # negatives
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccccc1", "benzene", "苯"),
    ("Cc1ccc2ccccc2c1", "2-methylnaphthalene", "2-甲基萘"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_anthracene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
