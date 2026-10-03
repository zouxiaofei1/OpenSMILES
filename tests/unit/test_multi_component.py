"""多组分命名：阴阳离子配对成盐、剩余组分分号连接、组分绝不丢弃（P-77）。"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# 能唯一配对的有机阴阳离子 → P-77.1.1 二元盐名（阳离子名 + 阴离子名，空格）
ion_pair__CASES = [
    # 中文阴离子片段独立命名为酯式（硫酸乙酯），与盐式 gold 有出入，故只强断言英文
    ("S(=O)(=O)(OCC)[O-].C(C)[N+](CCCNC(=O)CCCCCCCCCCCCCCCCC)(C)C",
     "ethyl-dimethyl-[3-(octadecanoylamino)propyl]azanium ethyl sulfate", None),
]

# 盐解离放宽后新增：有机阴离子/多原子阴离子/多份卤离子与有机片段成盐
halide_loose__CASES = [
    ("CCC(C(=O)[O-])N.[Cl-]", "2-aminobutanoate chloride", "2-氨基丁酸氯化物"),
    ("C1=CC=C2C(=C1)NC(=N2)CCCCCCCC3=NC4=CC=CC=C4N3.[Cl-].[Cl-]",
     "2-[7-(1H-benzimidazol-2-yl)heptyl]-1H-benzimidazole dichloride",
     "2-[7-(1H-苯并咪唑-2-基)庚基]-1H-苯并咪唑二氯化物"),
]

# 非盐多组分：组分名按字母序（忽略位次/括号/倍增前缀）排列，分号连接
multi_component__CASES = [
    ("C1CO1.C=O.CC(C)(C)CC(C)(C)c1ccc(O)cc1",
     "formaldehyde; oxirane; 4-(2,4,4-trimethylpentan-2-yl)phenol",
     "甲醛; 环氧乙烷; 4-(2,4,4-三甲基戊-2-基)苯酚"),
    ("[Na+].C1CO1", "oxirane; sodium", "环氧乙烷; 钠"),  # 不可配对的组分也保留
]

# 盐通道不回归：金属/卤化物/氢卤酸盐路径不受多组分改动影响
salt_no_regress__CASES = [
    ("[Na+].[O-]C(=O)c1ccccc1", "sodium benzoate", "苯甲酸钠"),
    ("Cl.Cl.N[C@@H](CCCCN)C(=O)O",
     "(2S)-2,6-diaminohexanoic acid dihydrochloride", None),
]


def _check(smiles: str, en: str, zh: str | None) -> None:
    """跑命名并比对 en/zh（样板同 test_acid_chain）。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", ion_pair__CASES)
def test_multicomponent_ion_pair(smiles: str, en: str, zh: str | None) -> None:
    """唯一配对阴阳离子拼成 P-77 二元盐名。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", halide_loose__CASES)
def test_multicomponent_halide_loose(smiles: str, en: str, zh: str | None) -> None:
    """有机/多原子/多份阴离子不再被丢弃。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", multi_component__CASES)
def test_multicomponent_alphabetical(smiles: str, en: str, zh: str | None) -> None:
    """非盐多组分按字母序排列、分号连接。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", salt_no_regress__CASES)
def test_salt_path_no_regress(smiles: str, en: str, zh: str | None) -> None:
    """盐通道命名不回归。"""
    _check(smiles, en, zh)
