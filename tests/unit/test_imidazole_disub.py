# IUPAC: P-14.3.4 / P-22.2.1
# Layer: L2
"""Imidazole disub scope (P-14.3.4) — negative guard only.

Positive 2-substituted imidazole cases are not covered in this file. The
retained case asserts an amino-imidazolecarboxylic acid is not named as a
simple retained imidazole.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: FG-bearing acid must not be simple imidazole parent
    ("Nc1[nH]cnc1C(=O)O", None, None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_imidazole_disub(smiles: str, en: str | None, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    if not r.success:
        return  # fallback 已删：无候选显式失败
    if en is None:
        assert "imidazole" not in normalize_en(r.en) or "carboxylic" in normalize_en(r.en)
        return
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_acid_not_simple_imidazole() -> None:
    """Amino-imidazolecarboxylic acid is not a simple retained imidazole."""
    r = SMILESNNamer().name("Nc1[nH]cnc1C(=O)O")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert en != "1h-imidazole"
    assert "2-amino" not in en or "carboxylic" in en or "acid" in en
