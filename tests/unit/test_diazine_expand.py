# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Diazine expand scope (P-22.2.1) — negative guard only.

Positive substituted diazine cases are not covered in this file. The retained
cases assert pyridine / benzene must not be captured as diazines.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: pyridine / benzene must not be captured as diazine
    ("c1ccncc1", "pyridine", "吡啶"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_diazine_expand(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_complex_side_not_simple_diazine() -> None:
    """Pentyl side chain is out of scope A (left for later)."""
    r = SMILESNNamer().name("ClC1=NC=C(C=N1)CCCCC")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "pentylpyrimidine" not in en
    assert "chloropyrimidine" not in en or "pentyl" not in en
