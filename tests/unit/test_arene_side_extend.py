# IUPAC: P-29.3 / P-14.3.4
# Layer: L2,L3,L5
"""Arene simple side-chain extensions: ω-haloalkyl, alkoxy chains, quinoline CF3.

P-29.3 / P-14.3.4: ring–(CH2)n–X ω-halo n-alkyl as compound substituent.
P-22.2.1: quinoline with halo + CF3; extended alkoxy (2-methoxyethoxy).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # A: ω-halo n-alkyl on benzene
    ("BrCCCC1=CC=CC=C1", "(3-bromopropyl)benzene", "(3-溴丙基)苯"),
    # negatives: keep existing correct behaviour
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_side_extend(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
