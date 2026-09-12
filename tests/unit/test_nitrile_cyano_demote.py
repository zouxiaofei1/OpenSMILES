# IUPAC: P-61.1.3,P-66.1.1
# Layer: L1,L2,L3
"""非主官能团的腈须前缀化为 cyano，不得被词干链吞碳把 N 悬空误命名成 amino。

腈碳(C≡N)是「自带碳的前缀叶」：当自由基(p41=1)/羧酸/酰胺等更高优先级主基团存在时，
腈退出主基团并保留 cyano 叶身份——其碳须从开链主链剔除（同中性 COOH 的 carboxy 叶，
P-61.1.3），否则如 `N#CCC*`(→*CH₂CH₂CN) 会被读成饱和丙基链 + 端 N 而错名 3-aminopropyl；
`HOOC-CH₂CN` 也被错名成 3-aminopropanoic acid。

- 链端/支链腈同规则：`*CH₂CH₂CN` → 2-cyanoethyl；`HOOC-CH₂CN` → 2-cyanoacetic acid；
  `HOOC-CH(CN)CH₃` → 2-cyanopropanoic acid（腈碳不得占主链位）。
- 无更高优先级基团时腈仍是主基团后缀（butanenitrile/…），不得 cyano 化。
"""

from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer

# 自由基/前缀残基正例：端基 C≡N 不并入 -yl 词干。
POS_FRAG = [
    ("N#CCC*", "2-cyanoethyl", "2-氰基乙基"),
    ("N#CC*", "cyanomethyl", "氰基甲基"),
    ("N#CCCC*", "3-cyanopropyl", "3-氰基丙基"),
]


@pytest.mark.parametrize("smiles,en,zh", POS_FRAG)
def test_cyano_radical_fragment(smiles: str, en: str, zh: str) -> None:
    """锚定自由基残基：腈作前缀叶，N 不得误命名成 amino。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh


# 整分子正例：酸/酰胺主基团压制腈（链端或支链），腈碳不进主链。
POS_WHOLE = [
    ("OC(=O)CC#N", "2-cyanoacetic acid", "2-氰基乙酸"),
    ("OC(=O)CCC#N", "3-cyanopropanoic acid", "3-氰基丙酸"),
    ("CC(C#N)C(=O)O", "2-cyanopropanoic acid", "2-氰基丙酸"),
    ("OC(=O)C(C#N)CC", "2-cyanobutanoic acid", "2-氰基丁酸"),
    ("NC(=O)CC#N", "2-cyanoacetamide", "2-氰基乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", POS_WHOLE)
def test_cyano_whole_molecule(smiles: str, en: str, zh: str) -> None:
    """整分子：腈被更高优先级基团压制 → cyano 前缀叶，主基团后缀不变。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh


# 负例：无更高优先级基团时腈保持主基团（后缀 nitrile），不得 cyano 化 / amino 化。
NEG = [
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("NCC#N", "aminoacetonitrile", "氨基乙腈"),
    ("O=CCC#N", "3-oxopropanenitrile", "3-氧代丙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", NEG)
def test_nitrile_principal_untouched(smiles: str, en: str, zh: str) -> None:
    """腈作主基团（后缀 -nitrile/-腈）时不受 cyano 前缀化影响。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh
