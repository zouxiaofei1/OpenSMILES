# IUPAC: P-62.2 / P-31.1
# Layer: L2,L4,L5
"""Unsaturated primary amines (alkenamines).

Amines use the segment-style ene insertion (same as alcohols):
but-3-en-1-amine. Previously the amine entry had no ene segment, so double
bonds were silently dropped (butan-1-amine).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("C=CCCN", "but-3-en-1-amine", "丁-3-烯-1-胺"),
    ("CC=CCN", "but-2-en-1-amine", "丁-2-烯-1-胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenamine(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
