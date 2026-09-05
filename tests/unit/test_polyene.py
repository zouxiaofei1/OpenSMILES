# IUPAC: P-31.1
# Layer: L2,L4,L5
"""Unsubstituted acyclic polyenes (alkadiene / alkatriene).

Parent chain through all non-aromatic C=C; ene locant set lowest.
English: buta-1,3-diene; Chinese: 丁-1,3-二烯.
Monoalkene / alkyne negatives must stay unchanged.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic polyenes (≥2 C=C, no ring/FG/triple)
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
    ("C=CC=CC=C", "hexa-1,3,5-triene", "己-1,3,5-三烯"),
    ("C=CCC=C", "penta-1,4-diene", "戊-1,4-二烯"),
    ("C=CCCC=C", "hexa-1,5-diene", "己-1,5-二烯"),
    ("C=CC=CC", "penta-1,3-diene", "戊-1,3-二烯"),
    ("CC=CC=CC", "hexa-2,4-diene", "己-2,4-二烯"),
    ("C=CCCCC=C", "hepta-1,6-diene", "庚-1,6-二烯"),
    # negative: monoalkene / alkyne must not become polyene
    ("CC=CC", "but-2-ene", "丁-2-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
