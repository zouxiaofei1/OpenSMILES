# IUPAC: P-63.2.1/P-63.2.2  (醇→烷氧基、硫醇→烷基硫基；简单取代基+氧/硫连接无需围栏)
# Layer: L3
"""简单取代基 + 氧/硫 (R-O-/R-S- 取代基前缀) 不加括号。

- 规则：当取代基是"简单取代基 + -oxy/-sulfanyl 桥"（如 propan-2-yloxy、
  hexadecanoyloxy、acetyloxy、benzyloxy、naphthalen-1-yloxy、propan-2-ylsulfanyl），
  且前端(R)本身是不需括号的简单取代基时，整个前缀不加围栏（ChEBI/gold 平铺式）。
- 反例：前端被取代（4-nitrophenyl-、2,6-dichlorophenyl-）时前缀须加围栏；括号只括前端、
  -oxy 留括号外（5-(trifluoromethyl)pyridin-2-yl]oxyphenoxy，P-63.2.2.1.1 的 (pyridin-2-yl)oxy）。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# (smiles, expected_en, expected_zh_or_None) —— EN 期望即 gold(ChEBI) 无括号写法
CASES = [
    # acyloxy：前端 hexadecanoyl / acetyl 简单 → 不加括号
    (
        "CCCCCCCCCCCCCCCC(=O)OC(CCCCCCCCC)CCCCCCCC(=O)O",
        "9-hexadecanoyloxyoctadecanoic acid",
        None,
    ),
    (
        "CCCCCCCCCCCCCCCc1ccc(C(=O)O)c(OC(C)=O)c1",
        "2-acetyloxy-4-pentadecylbenzoic acid",
        None,
    ),
    # alkoxy：前端 propan-2-yl 简单（带自由价位次但仍简单）→ 不加括号
    (
        "BrC=1C=C(C=O)C=C(C1)OC(C)C",
        "3-bromo-5-propan-2-yloxybenzaldehyde",
        None,
    ),
    (
        "C(Oc1ccccc1)c1ccccc1",
        "benzyloxybenzene",
        "苄氧基苯",
    ),
    # sulfanyl：前端 propan-2-yl 简单 → 不加括号
    (
        "CC(C)SCC#N",
        "2-propan-2-ylsulfanylacetonitrile",
        None,
    ),
    # 糖苷/手性环基前端 naphthalen-1-yl 简单，thiophen-2-yl 亦不括（gold 写法）
    (
        "CNCC[C@H](Oc1cccc2ccccc12)c1cccs1.Cl",
        "(3S)-N-methyl-3-naphthalen-1-yloxy-3-thiophen-2-ylpropan-1-amine hydrochloride",
        None,
    ),
    # sulfooxy：前端 sulfo 为 retained 简单前缀 → 不加括号
    (
        "N[C@@H](Cc1ccc(OS(=O)(=O)O)cc1)C(=O)O",
        "(2S)-2-amino-3-(4-sulfooxyphenyl)propanoic acid",
        None,
    ),
]

# 负例：复杂前端或非氧/硫连接，括号必须保留（原样不回退）
NEGATIVES = [
    # 前端 4-nitrophenyl 被取代 → 整体加括号
    ("[N+](=O)([O-])C1=CC=C(C=C1)SCCO", "2-(4-nitrophenylsulfanyl)ethanol"),
    # 前端 2,6-dichlorophenyl 被取代 → 整体加括号
    (
        "Clc1ccc(CCC(Cn2ccnc2)Sc2c(Cl)cccc2Cl)cc1",
        "1-[4-(4-chlorophenyl)-2-(2,6-dichlorophenylsulfanyl)butyl]imidazole",
    ),
    # 前端 5-(trifluoromethyl)pyridin-2-yl 被取代 → 前端括起、oxy 留括号外（P-63.2.2.1.1：
    # (pyridin-2-yl)oxy；P-65.6.3.2.3：3-[(pyridine-3-carbonyl)oxy]…）
    (
        "CCCCOC(=O)C(C)Oc1ccc(Oc2ccc(C(F)(F)F)cn2)cc1",
        "butyl 2-[4-[5-(trifluoromethyl)pyridin-2-yl]oxyphenoxy]propanoate",
    ),
    # -amino 连接不在氧/硫规则内：benzylamino 保留括号
    ("OCCNCc1ccccc1", "2-(benzylamino)ethanol"),
    # retained 简单烷氧基本就无括号（methoxy/ethoxy），不能回归成加括号
    ("CCOc1ccccc1", "ethoxybenzene"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_obridge_simple_no_paren(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en", NEGATIVES)
def test_obridge_complex_keeps_paren(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
