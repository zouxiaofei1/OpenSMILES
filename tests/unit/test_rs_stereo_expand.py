# IUPAC: P-92 / P-93
# Layer: L5
"""R/S coverage expand: ketone/ester/sat-hetero + ester slot + collapse gate."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # ester: RS between alkyl and acyl (not before methyl)
    (
        "COC(=O)[C@@H](N)CC(C)C",
        "methyl (2S)-2-amino-4-methylpentanoate",
        "(2S)-2-氨基-4-甲基戊酸甲酯",
    ),
    # ketone multi-center
    (
        "O=C(CO)[C@@H](O)[C@H](O)CO",
        "(3S,4R)-1,3,4,5-tetrahydroxypentan-2-one",
        "(3S,4R)-1,3,4,5-四羟基戊-2-酮",
    ),
    # ketone single center with locant
    (
        "CC(=O)[C@H](O)c1ccccc1",
        "(1R)-1-hydroxy-1-phenylpropan-2-one",
        "(1R)-1-羟基-1-苯基丙-2-酮",
    ),
    # regression open-chain acid
    (
        "O=C(O)[C@@H](O)CO",
        "(2S)-2,3-dihydroxypropanoic acid",
        "(2S)-2,3-二羟基丙酸",
    ),
    # regression E/Z + R/S merge
    (
        "C[C@@H](O)CCCC/C=C/C(=O)O",
        "(2E,8R)-8-hydroxynon-2-enoic acid",
        "(2E,8R)-8-羟基壬-2-烯酸",
    ),
    # negatives: no spurious R/S
    ("CCCCCCCCCCCC(=O)O", "dodecanoic acid", "十二酸"),
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_rs_stereo_expand(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
