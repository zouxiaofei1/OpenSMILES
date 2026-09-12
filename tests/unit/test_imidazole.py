# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained 1H-imidazole scope (P-22.2.1) — negative guard only.

Positive imidazole cases are not covered in this file. The retained cases
assert benzene / pyridine / open chain must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: must not regress mono-hetero5 / arene / open chain
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CCCCC", "pentane", "戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_imidazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
