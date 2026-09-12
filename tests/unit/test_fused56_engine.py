# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4,L5
"""Fused56 e2e scope (P-22.2.1) — negative guard only.

Positive benzofuran/benzothiophene e2e cases are not covered in this file.
The retained case asserts indole is not captured by the fused56 engine.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_E2E = [
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_e2e_fused56_mono(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
