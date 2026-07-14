# IUPAC: P-29.3 / P-14.3.4
# Layer: L2,L3,L5
"""Arene simple side-chain extensions: ω-haloalkyl, alkoxy chains, quinoline CF3.

P-29.3 / P-14.3.4: ring–(CH2)n–X ω-halo n-alkyl as compound substituent.
P-22.2.1: quinoline with halo + CF3; extended alkoxy (2-methoxyethoxy).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # A: ω-halo n-alkyl on benzene
    ("BrCCCC1=CC=CC=C1", "(3-bromopropyl)benzene", "(3-溴丙基)苯"),
    # B: extended alkoxy (2-methoxyethoxy) + ring halo
    (
        "BrC1=CC(=C(C=C1)OCCOC)F",
        "4-bromo-2-fluoro-1-(2-methoxyethoxy)benzene",
        "4-溴-2-氟-1-(2-甲氧基乙氧基)苯",
    ),
    # C: quinoline halo + CF3
    (
        "BrC1=CC=C2C=CC(=NC2=C1)C(F)(F)F",
        "7-bromo-2-(trifluoromethyl)quinoline",
        "7-溴-2-三氟甲基喹啉",
    ),
    (
        "ClC1=NC2=CC=CC=C2C(=C1)C(F)(F)F",
        "2-Chloro-4-(trifluoromethyl)quinoline",
        "2-氯-4-三氟甲基喹啉",
    ),
    # negatives: keep existing correct behaviour
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    ("FC(F)(F)c1ccccc1", "trifluoromethylbenzene", None),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_side_extend(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
