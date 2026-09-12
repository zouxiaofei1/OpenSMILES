# IUPAC: P-62.2.1 / P-63.1.4 / P-22.2.1
# Layer: L2,L4,L5
"""Pyridin-amine / pyridin-ol scope (P-62.2.1) — negative guard only.

Positive pyridinamine/pyridinol cases are not covered in this file. The retained
cases assert bare pyridine, aniline, phenol, chain amine stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negative: bare pyridine, aniline, phenol, chain amine, pyridinecarboxylic
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CCN", "ethanamine", "乙胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinamine_ol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
