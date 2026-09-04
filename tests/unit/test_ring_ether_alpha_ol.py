# IUPAC: P-63.1.1 / P-22.1.1
# Layer: L2 (parent_ownership) + L5
"""环醚 α-羟基(半缩醛型)醇:醇碳同时邻接另一单键氧(环醚/醚 O)时,醇 OH 氧不得游离出
owned_atoms 被 coverage 补齐阶段命名成假 hydroxy 前缀(与 -ol 后缀双算)。根因:_single_o_idx
盲取醇碳首个单键 O,歧义时取到醚氧。"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh") — 醇碳邻另一单键 O 的饱和杂环/醚醇
CASES = [
    # 核心 bug:THF-2-醇(oxolan-2-ol),2 位碳邻环醚 O + 羟基 O 两个单键氧
    (
        "C1OC(O)CC1",
        "oxolan-2-ol",
        "四氢呋喃-2-醇",
    ),
    # 同场景 6 元环:oxan-2-ol(THP-2-醇)
    (
        "C1OC(O)CCC1",
        "oxan-2-ol",
        "氧杂环己烷-2-醇",
    ),
    # acyclic 醚-醇(半缩醛/缩醛型):CH3-O-CH2-OH,醇碳邻醚 O + 羟基 O
    (
        "COCO",
        "methoxymethanol",
        "甲氧基甲醇",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ether_alpha_ol_no_bogus_hydroxy(smiles: str, en: str, zh: str) -> None:
    """醇氧属母体,不得再命名出 (hydroxy) 前缀。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# negatives — 不受该歧义影响、不应回归的相邻场景
@pytest.mark.parametrize(
    "smiles,en,zh",
    [
        ("C1CCC(O)CC1", "cyclohexanol", "环己醇"),
        ("OC1CCCC1", "cyclopentanol", "环戊醇"),
        ("OCC(O)C", "propane-1,2-diol", "丙烷-1,2-二醇"),
    ],
)
def test_unaffected_alcohols_still_fine(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
