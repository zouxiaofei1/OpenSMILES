# IUPAC: P-91.2 / P-93.4
# Layer: L5
"""E/Z stereo prefixes for alkenal, alkenenitrile, alkene, and polyene parents."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positives — mono-ene FG with BondStereo
    (
        r"CCCCCC/C=C/C=O",
        "(E)-non-2-enal",
        "(E)-壬-2-烯醛",
    ),
    (
        r"CCCCCCC/C=C/C=O",
        "(E)-dec-2-enal",
        "(E)-癸-2-烯醛",
    ),
    (
        r"CCCCCCCC/C=C\CCCCCCCC#N",
        "(Z)-octadec-9-enenitrile",
        "(Z)-十八-9-烯腈",
    ),
    # alkene parent
    (
        r"C/C=C/C",
        "(E)-but-2-ene",
        "(E)-丁-2-烯",
    ),
    (
        r"C/C=C\C",
        "(Z)-but-2-ene",
        "(Z)-丁-2-烯",
    ),
    # polyene: terminal C=C has no E/Z; only internal stereo should prefix
    (
        r"C=CCCCCCCCCCCC/C=C\CCCC",
        "(14Z)-nonadeca-1,14-diene",
        "(14Z)-十九-1,14-二烯",
    ),
    # negatives — no forced stereo → no (E)/(Z)
    ("C=CCC=O", "but-3-enal", "丁-3-烯醛"),
    ("C=CCC", "but-1-ene", "丁-1-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ez_stereo(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
