# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained 1H-pyrazole scope (P-22.2.1) — negative guard only.

Positive pyrazole cases are not covered in this file. The retained case asserts
pyridine is not misclassified as pyrazole.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: must not misclassify 1,3-diazole / diazine / mono-hetero5
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyrazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
