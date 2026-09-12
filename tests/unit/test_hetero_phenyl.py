# IUPAC: P-22.2.1 / P-29.3 / P-14.3.4
# Layer: L2,L3,L5
"""Phenyl-on-heteroarene scope (P-22.2.1) — negative guard only.

Positive phenyl-diazine / phenyl-hetero5 cases are not covered in this file.
The retained case asserts pyridine is not captured as a phenyl-hetero parent.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # regressions
    ("c1ccncc1", "pyridine", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hetero_phenyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
