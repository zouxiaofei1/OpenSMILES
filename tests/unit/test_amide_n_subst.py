# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_n_alkyl_amide.py: Open-chain N-substituted amide scope (P-66.1.1.1.3) — negative guard only.
test_amide_n_benzyl.py: Amide N-benzyl (and substituted benzyl) claim — not N-methyl.
test_n_phenyl_amide.py: Open-chain N-phenyl amide scope (P-66.1.1.1.3) — negative guard only.
test_n_hydroxy_amide.py: N-羟基(异羟肟酸型)酰胺 scope — 正例 + 近邻负例。
test_n_heteroaryl_benzamide.py: N-heteroaryl benzamide via recursive substituent path.
test_n_alkyl_amino_alcohol.py: N-烷基取代氨基醇的递归命名：仲胺带链外烷基侧链时不认领为
test_recursive_benzamide_north_star.py:
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_n_alkyl_amide.py
# IUPAC: P-66.1.1.1.3
# Layer: L1,L2,L3,L5
#
# Open-chain N-substituted amide scope (P-66.1.1.1.3) — negative guard only.
#
# Positive N-alkyl amide cases are not covered in this file. The retained cases
# assert primary amide / acid / aldehyde are not named as N-alkyl amides.
# ==========================================================================
n_alkyl_amide__CASES = [
    ("CC=O", "acetaldehyde", "乙醛"),
]


@pytest.mark.parametrize("smiles,en,zh", n_alkyl_amide__CASES)
def test_n_alkyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_amide_n_benzyl.py
# IUPAC: P-66.1.1.1.3 / P-29.3
# Layer: L2,L3
#
# Amide N-benzyl (and substituted benzyl) claim — not N-methyl.
#
# Open-chain monoamide parent; N–CH2–Ph (Ph may carry simple leaves) is
# N-benzyl / N-(…benzyl). N-phenyl and plain N-alkyl stay correct.
# ==========================================================================
amide_n_benzyl__CASES = [
    # negatives: N-phenyl, plain N-alkyl, primary amide
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", amide_n_benzyl__CASES)
def test_amide_n_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_n_phenyl_amide.py
# IUPAC: P-66.1.1.1.3
# Layer: L2,L3,L5
#
# Open-chain N-phenyl amide scope (P-66.1.1.1.3) — negative guard only.
#
# Positive N-phenyl amide cases are not covered in this file. The retained cases
# assert primary amide and C-phenyl benzamide are not named as N-phenyl amides.
# ==========================================================================
n_phenyl_amide__CASES = [
    ("c1ccc(C(=O)N)cc1", "benzamide", "苯甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", n_phenyl_amide__CASES)
def test_n_phenyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_n_hydroxy_amide.py
# IUPAC: P-66.1.1.1.3
# Layer: L1,L2,L3,L5
#
# N-羟基(异羟肟酸型)酰胺 scope — 正例 + 近邻负例。
#
# N-羟基酰胺/异羟肟酸: 羰基 C(=O) 上的酰胺 N 除羰基碳外还带一个不再连碳的
# 羟基氧(N-OH),如 acetohydroxamic acid 按 N-hydroxy-…amide 表达。此前 L1
# `_amide_n_info` 只允许酰胺 N 上的取代基为 C/H,N-OH 使酰胺漏检、羰基碳被误判
# 成醛 → 输出 1-[…]acetaldehyde(P-66.1.1.1.3 N-前缀命名)。
#
# 负例: N-烷基酰胺、伯酰胺、酮、酸、醛不得被误伤。
# ==========================================================================
n_hydroxy_amide__CASES = [
    # positive: N-hydroxy amide family (gold chebi-1411)
    ("CC(=O)N(O)CCCN", "N-(3-aminopropyl)-N-hydroxyacetamide", "N-(3-氨基丙基)-N-羟基乙酰胺"),
    ("CC(=O)NO", "N-hydroxyacetamide", "N-羟基乙酰胺"),
    ("CC(=O)N(O)C", "N-hydroxy-N-methylacetamide", "N-羟基-N-甲基乙酰胺"),
    ("CC(=O)N(C)O", "N-hydroxy-N-methylacetamide", "N-羟基-N-甲基乙酰胺"),
    ("CC(=O)N(O)O", "N,N-dihydroxyacetamide", "N,N-二羟基乙酰胺"),
    # negative: N-alkyl amide / primary amide / ketone / acid / aldehyde stay correct
    ("CC(=O)NCCCN", "N-(3-aminopropyl)acetamide", "N-(3-氨基丙基)乙酰胺"),
    ("CC(=O)NC", "N-methylacetamide", "N-甲基乙酰胺"),
    ("CC(=O)CCC", "pentan-2-one", "戊-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", n_hydroxy_amide__CASES)
def test_n_hydroxy_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_n_heteroaryl_benzamide.py
# IUPAC: P-66.1.1 / P-29 / P-22.2.1
# Layer: L2,L3,L5
#
# N-heteroaryl benzamide via recursive substituent path.
# ==========================================================================
n_heteroaryl_benzamide__REGRESS = [
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", n_heteroaryl_benzamide__REGRESS)
def test_n_heteroaryl_benzamide_regress(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_n_alkyl_amino_alcohol.py
# IUPAC: P-65.1 / P-29.3；Layer: L2,L3,L5
#
# N-烷基取代氨基醇的递归命名：仲胺带链外烷基侧链时不认领为
# "amino"，整条侧链交由递归后端命名（防截断成 aminoethanol）。
# ==========================================================================
n_alkyl_amino_alcohol__CASES = [
    # 断点 case：仲胺 N 同时连链碳与侧链碳（深度递归两层）
    (
        "OCCNCCNCCC",
        "2-[2-(propylamino)ethylamino]ethanol",
        "2-[2-(丙氨基)乙氨基]乙醇",
    ),
    # 同一分子反转输入：结果必须一致（此前依赖输入原子序）
    (
        "CCCNCCNCCO",
        "2-[2-(propylamino)ethylamino]ethanol",
        "2-[2-(丙氨基)乙氨基]乙醇",
    ),
    # N-甲基 / N-氨基乙基（自由基链 1 位省略：2-aminoethylamino）
    ("OCCNCCN", "2-(2-aminoethylamino)ethanol", "2-(2-氨基乙氨基)乙醇"),
    # N-苄基
    ("OCCNCc1ccccc1", "2-(benzylamino)ethanol", "2-(苄氨基)乙醇"),
    # N 桥接二醇（母体链不含 N，两侧都是侧链）
    ("OCCNCCO", "2-(2-hydroxyethylamino)ethanol", "2-(2-羟基乙氨基)乙醇"),

    # 大分子保持多层递归（饱和链自由价 1 位省略：ethylamino/pentylamino）
    (
        "NCCNCCCCCNCCNCCNCCCNCCCNCCO",
        "2-[3-[3-[2-[2-[5-(2-aminoethylamino)pentylamino]"
        "ethylamino]ethylamino]propylamino]"
        "propylamino]ethanol",
        "2-[3-[3-[2-[2-[5-(2-氨基乙氨基)戊氨基]"
        "乙氨基]乙氨基]丙氨基]丙氨基]乙醇",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", n_alkyl_amino_alcohol__CASES)
def test_n_alkyl_amino_alcohol(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_recursive_benzamide_north_star.py
# tests/unit/test_recursive_benzamide_north_star.py
# IUPAC: P-66.1 / P-29 / P-25 / P-65
# Layer: L2–L5
# ==========================================================================
recursive_benzamide_north_star__NORTH = (
    "COc1ccc(C(=O)Nc2nc3c(s2)CCC3C(=O)N2CCCC2)cc1",
    "4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide",
    "4-甲氧基-N-(4-(吡咯烷-1-羰基)-5,6-二氢-4H-环戊并[d]噻唑-2-基)苯甲酰胺",
)
