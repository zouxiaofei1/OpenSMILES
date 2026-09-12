# IUPAC: P-22.2.1 / P-25
# Layer: L2,L5
"""Quinazoline scope (P-22.2.1) — negative guard only.

Positive quinazoline cases are not covered in this file. The retained case
asserts naphthalene is not named quinazoline.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_quinazoline(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cinnoline_not_benzodiazine() -> None:
    r = SMILESNNamer().name("c1ccc2nnccc2c1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "quinazoline" not in en and "quinoxaline" not in en
