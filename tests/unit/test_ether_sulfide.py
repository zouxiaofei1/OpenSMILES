# 合并自 5 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_dialkyl_ether.py: Open-chain simple dialkyl ethers (symmetric retained + alkoxyalkane).
test_branched_ether.py: Branched dialkyl ethers: isopropyl / hexafluoroisopropyl functional class.
test_dialkyl_sulfide.py: Open-chain dialkyl sulfide scope (P-63.2.1) — negative guard only.
test_anisole.py: Aromatic alkoxy prefixes (methoxy/ethoxy) and retained anisole (P-63.2.2).
test_obridge_simple_no_paren.py: 简单取代基 + 氧/硫 (R-O-/R-S- 取代基前缀) 不加括号。
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_dialkyl_ether.py
# IUPAC: P-63.2.2
# Layer: L1,L2,L5
#
# Open-chain simple dialkyl ethers (symmetric retained + alkoxyalkane).
#
# Exactly one ether O (two C neighbors, no H, not ester/anhydride). Both arms
# unsubstituted linear alkyl C1–C4. Symmetric: diethyl ether / 二乙基醚.
# Asymmetric: methoxyethane / 甲氧基乙烷; 1-methoxypropane / 1-甲氧基丙烷.
# ==========================================================================
dialkyl_ether__CASES = [
    # positive: asymmetric alkoxyalkane
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCCOC", "1-methoxypropane", "1-甲氧基丙烷"),
    ("CCOCCC", "1-ethoxypropane", "1-乙氧基丙烷"),
    ("C=CC(=O)OC", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
]


@pytest.mark.parametrize("smiles,en,zh", dialkyl_ether__CASES)
def test_dialkyl_ether(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_branched_ether.py
# IUPAC: P-63.2.2
# Layer: L2, L5
#
# Branched dialkyl ethers: isopropyl / hexafluoroisopropyl functional class.
#
# Linear C1–C4 arms keep alkoxyalkane / retained sym ether. Branched specials
# use functional-class names (… methyl ether / …基甲醚).
# ==========================================================================
branched_ether__CASES = [
    # hexafluoroisopropyl methyl ether (benchmark tiers-15987)
    (
        "COC(C(F)(F)F)C(F)(F)F",
        "hexafluoroisopropyl methyl ether",
        "六氟异丙基甲醚",
    ),
    # unsubstituted isopropyl methyl ether
    ("COC(C)C", "isopropyl methyl ether", "异丙基甲醚"),
    # linear regressions
    ("COC", "dimethyl ether", "二甲基醚"),
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCOCC", "diethyl ether", "二乙基醚"),
    # negatives
    ("CCO", "ethanol", "乙醇"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
]


# ==========================================================================
# 合并自 test_dialkyl_sulfide.py
# IUPAC: P-63.2.1
# Layer: L1,L2,L5
#
# Open-chain dialkyl sulfide scope (P-63.2.1) — negative guard only.
#
# Positive sulfide cases are not covered in this file. The retained cases
# assert ethanol / sulfoxide / thioester are not named as dialkyl sulfides.
# ==========================================================================
dialkyl_sulfide__CASES = [
    # negative: thiol / ether / alcohol must not break
    ("CCO", "ethanol", "乙醇"),
]

# sulfoxide / thioester must not be named as dialkyl sulfide
dialkyl_sulfide__NEG_NOT_SULFIDE = [
    ("CS(C)=O", "dimethyl sulfide"),
    ("CC(=O)SC", "methyl sulfide"),
]


@pytest.mark.parametrize("smiles,forbidden_en", dialkyl_sulfide__NEG_NOT_SULFIDE)
def test_not_dialkyl_sulfide(smiles: str, forbidden_en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert en != normalize_en(forbidden_en)
    assert "sulfide" not in en
    assert "硫醚" not in normalize_zh(r.zh)


@pytest.mark.parametrize("smiles,en,zh", dialkyl_sulfide__CASES)
def test_dialkyl_sulfide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_anisole.py
# IUPAC: P-63.2.2
# Layer: L2,L3,L5
#
# Aromatic alkoxy prefixes (methoxy/ethoxy) and retained anisole (P-63.2.2).
#
# Parent = benzene or phenol; ring alkoxy O–R (R = Me/Et linear) is a prefix.
# Unsubstituted methoxybenzene → retained anisole / 甲氧基苯（对齐 benchmark 金标）.
# Do not absorb open-chain ethers or bare phenol/benzene.
# ==========================================================================
anisole__CASES = [
    # positive: methoxyphenol
    ("COc1ccc(O)cc1", "4-methoxyphenol", "4-甲氧基苯酚"),
    # positive: methoxy + halo / methyl (lowest set of locants)
    ("COc1ccccc1I", "1-iodo-2-methoxybenzene", None),
    ("COc1ccccc1C", "1-methoxy-2-methylbenzene", None),
]


@pytest.mark.parametrize("smiles,en,zh", anisole__CASES)
def test_anisole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_obridge_simple_no_paren.py
# IUPAC: P-63.2.1/P-63.2.2  (醇→烷氧基、硫醇→烷基硫基；简单取代基+氧/硫连接无需围栏)
# Layer: L3
#
# 简单取代基 + 氧/硫 (R-O-/R-S- 取代基前缀) 不加括号。
#
# - 规则：当取代基是"简单取代基 + -oxy/-sulfanyl 桥"（如 propan-2-yloxy、
#   hexadecanoyloxy、acetyloxy、benzyloxy、naphthalen-1-yloxy、propan-2-ylsulfanyl），
#   且前端(R)本身是不需括号的简单取代基时，整个前缀不加围栏（ChEBI/gold 平铺式）。
# - 反例：前端被取代（4-nitrophenyl-、2,6-dichlorophenyl-）时前缀须加围栏；括号只括前端、
#   -oxy 留括号外（5-(trifluoromethyl)pyridin-2-yl]oxyphenoxy，P-63.2.2.1.1 的 (pyridin-2-yl)oxy）。
# ==========================================================================
obridge_simple_no_paren__CASES = [
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
obridge_simple_no_paren__NEGATIVES = [
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


@pytest.mark.parametrize("smiles,en,zh", obridge_simple_no_paren__CASES)
def test_obridge_simple_no_paren(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en", obridge_simple_no_paren__NEGATIVES)
def test_obridge_complex_keeps_paren(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
