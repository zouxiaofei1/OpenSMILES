# IUPAC: P-92 / P-93
# Layer: L5
"""R/S CIP on expanded parent kinds: amide, nitrile, aldehyde, thiol, diacid."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positives — RDKit CIP calibrated
    (
        "CC[C@H](C)C(=O)N",
        "(2S)-2-methylbutanamide",
        "(2S)-2-甲基丁酰胺",
    ),
    (
        "CC[C@H](C)C#N",
        "(2S)-2-methylbutanenitrile",
        "(2S)-2-甲基丁腈",
    ),
    (
        "CC[C@H](C)C=O",
        "(2S)-2-methylbutanal",
        "(2S)-2-甲基丁醛",
    ),
    (
        "CC[C@H](C)S",
        "(2S)-butane-2-thiol",
        "(2S)-丁-2-硫醇",
    ),
    (
        "N[C@H](CCCC(=O)O)C(=O)O",
        "(2R)-2-aminohexanedioic acid",
        "(2R)-2-氨基己二酸",
    ),
    (
        "FC1=CC=C(C=C1)[C@@H](CC(=O)N)C=C",
        "(3S)-3-(4-fluorophenyl)pent-4-enamide",
        "(3S)-3-(4-氟苯基)戊-4-烯酰胺",
    ),
    # negatives — no spurious R/S
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCCCCCCCCCCC(=O)O", "dodecanoic acid", "十二酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_rs_stereo_kinds(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
