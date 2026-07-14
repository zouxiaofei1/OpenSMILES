# IUPAC: P-22.2.1 / P-25
# Layer: L2,L5
"""Unsubstituted quinazoline retained parent (benzodiazine 1,3)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("c1ccc2ncncc2c1", "quinazoline", "喹唑啉"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccc2ncccc2c1", "quinoline", "喹啉"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_quinazoline(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_quinoxaline_not_quinazoline() -> None:
    r = SMILESNNamer().name("c1ccc2nccnc2c1")
    assert r.success
    assert "quinazoline" not in normalize_en(r.en)
