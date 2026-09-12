# IUPAC: P-92 / P-93
# Layer: L4,L5
"""R/S CIP stereodescriptor prefixes on parent-chain chiral carbons."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positives: skeleton correct, only missing R/S
    (
        "O=C(O)[C@@H](O)CO",
        "(2S)-2,3-dihydroxypropanoic acid",
        "(2S)-2,3-二羟基丙酸",
    ),
    (
        "CC[C@H](C)C(=O)[O-]",
        "(2S)-2-methylbutanoate",
        "(2S)-2-甲基丁酸根",
    ),
    (
        "NCC[C@@H](O)C[C@H](N)C(=O)O",
        "(2S,4R)-2,6-diamino-4-hydroxyhexanoic acid",
        "(2S,4R)-2,6-二氨基-4-羟基己酸",
    ),
    (
        "O=C(O)C[C@H](O)CCCCCCCO",
        "(3R)-3,10-dihydroxydecanoic acid",
        "(3R)-3,10-二羟基癸酸",
    ),
    (
        "CCCCCCCCCCCC[C@H](O)C(=O)[O-]",
        "(2S)-2-hydroxytetradecanoate",
        "(2S)-2-羟基十四酸根",
    ),
    (
        "C[C@@H](O)CCCCCCC(=O)O",
        "(8R)-8-hydroxynonanoic acid",
        "(8R)-8-羟基壬酸",
    ),
    # E/Z + R/S merge (gold: (E,8R)-)
    (
        "C[C@@H](O)CCCC/C=C/C(=O)O",
        "(2E,8R)-8-hydroxynon-2-enoic acid",
        "(2E,8R)-8-羟基壬-2-烯酸",
    ),
    # 取代基/自由基母体：手性环成为取代基前缀时也应携带自身 R/S
    # （修复1：implicit-H 的 [C@]/[C@@] 被 RDKit 判为 3 配位而非立体中心）
    (
        "O[C@@]1[C@@](C*)OC(O)[C@](N)[C@]1O",
        "[(2R,3R,4S,5R)-5-amino-3,4,6-trihydroxyoxan-2-yl]methyl",
        "[(2R,3R,4S,5R)-5-氨基-3,4,6-三羟基氧杂环己烷-2-基]甲基",
    ),
    # negatives: no spurious R/S
    ("CCCCCCCCCCCC(=O)O", "dodecanoic acid", "十二酸"),
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    (
        "OC/C=C/C(=O)O",
        "(2E)-4-hydroxybut-2-enoic acid",
        "(2E)-4-羟基丁-2-烯酸",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_rs_stereo(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
