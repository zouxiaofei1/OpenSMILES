# IUPAC: P-41（HOOC-CH₂-CH₂- → 2-carboxyethyl 类比）/ P-61.1.3
# Layer: L1,L2,L3
"""端位 C≡N 在更高优先级主基团（自由基/羧酸/酰胺/酯…）下退为 cyano 前缀叶：

词干链不得吞掉腈碳——否则开链骨架把腈碳当饱和烷基碳、N 悬空误命名成 amino
（N#CCC* 曾错成 3-aminopropyl，应为 2-cyanoethyl）。修复点：L1
_arbitrate_parts 腈降级进 demoted_nitriles（同羧酸叶），L2 parent_skeleton
把降级腈碳剔除出开链主链。腈自身为主基团（无更高类）时不受影响。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# 片段级：* 锚定残基（L3 as_substituent 递归命名同一路径）。
POS_FRAG = [
    ("N#CC*", "cyanomethyl", "氰基甲基"),
    ("N#CCC*", "2-cyanoethyl", "2-氰基乙基"),
    ("N#CCCC*", "3-cyanopropyl", "3-氰基丙基"),
    ("N#CCCCC*", "4-cyanobutyl", "4-氰基丁基"),
]


@pytest.mark.parametrize("smiles,en,zh", POS_FRAG)
def test_cyano_alkyl_fragment(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子：腈作 cyano 前缀挂在酸/酰胺/酯主基上。
POS_WHOLE = [
    ("N#CCC(=O)O", "2-cyanoacetic acid", "2-氰基乙酸"),
    ("N#CCCC(=O)O", "3-cyanopropanoic acid", "3-氰基丙酸"),
    ("NC(=O)CC#N", "2-cyanoacetamide", "2-氰基乙酰胺"),
    ("N#CCC(=O)OC", "methyl cyanoacetate", "氰基乙酸甲酯"),
    ("N#CCCc1ccccc1C(=O)O", "2-(2-cyanoethyl)benzoic acid", "2-(2-氰基乙基)苯甲酸"),
    ("N#CCCc1ccc(C(=O)O)cc1", "4-(2-cyanoethyl)benzoic acid", "4-(2-氰基乙基)苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", POS_WHOLE)
def test_cyano_leaf_whole_molecule(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 近邻负例：腈仍是主基团（无更高类）或芳环 cyano 前缀，不得误伤。
NEG = [
    # 醇比腈低优先，腈保持主基。
    ("N#CCO", "hydroxyacetonitrile", "羟基乙腈"),
    # 醛（p41 高于腈）压成 oxo，腈保持主基。
    ("N#CCCCC=O", "5-oxopentanenitrile", "5-氧代戊腈"),
    # 芳环 cyano 前缀（酸主基）不变。
    ("N#Cc1ccc(C(=O)O)cc1", "4-cyanobenzoic acid", "4-氰基苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", NEG)
def test_cyano_not_demoted(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
