# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_alkene.py: Simple acyclic monoalkenes: ane→ene; lowest double-bond locant.
test_mono_alkyne.py: Simple acyclic monoalkynes: ane→yne; lowest triple-bond locant.
test_polyene.py: Unsubstituted acyclic polyenes (alkadiene / alkatriene).
test_enyne_polyyne.py: 混合烯炔与多炔命名：ene 前 yne 后、词干插 a、多重键集合编号优先。
test_alkyne_fgs.py: Open-chain mono alkyne FG parents: alkynol / alkynal / alkynenitrile / alkynamide.
test_ring_cumulene.py: 环内累积二烯 (ring cumulene) 编号：C=C=C 嵌入环内。
test_exocyclic_methylidene.py: Exocyclic double bond leaves: parent=CH2 must keep the double bond (methylidene).
test_long_chain_stems.py: Long-chain open parent stems C11–C35 (numerical multipliers + alkane stems).
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mono_alkene.py
# IUPAC: P-31.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic monoalkenes: ane→ene; lowest double-bond locant.
#
# C2/C3 omit locant (ethene/propene); C≥4 use stem-loc-ene / 位次-烯.
# Parent chain must include the unique non-aromatic C=C.
# ==========================================================================
mono_alkene__CASES = [
    ("CC=CCC", "pent-2-ene", "戊-2-烯"),
    ("C=CCCCC", "hex-1-ene", "己-1-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_alkene__CASES)
def test_mono_alkene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_mono_alkyne.py
# IUPAC: P-31.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic monoalkynes: ane→yne; lowest triple-bond locant.
#
# C2 retained acetylene/乙炔; C3 propyne/丙炔 (omit locant);
# C≥4 use stem-loc-yne / 位次-炔. Parent chain includes unique C≡C.
# ==========================================================================
mono_alkyne__CASES = [
    # positive: acyclic monoalkynes
    ("C#CC", "propyne", "丙炔"),
    ("CC#CC", "but-2-yne", "丁-2-炔"),
    ("CC#CCC", "pent-2-yne", "戊-2-炔"),
    ("C#CCCCC", "hex-1-yne", "己-1-炔"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_alkyne__CASES)
def test_mono_alkyne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_polyene.py
# IUPAC: P-31.1
# Layer: L2,L4,L5
#
# Unsubstituted acyclic polyenes (alkadiene / alkatriene).
#
# Parent chain through all non-aromatic C=C; ene locant set lowest.
# English: buta-1,3-diene; Chinese: 丁-1,3-二烯.
# Monoalkene / alkyne negatives must stay unchanged.
# ==========================================================================
polyene__CASES = [
    ("C=CC=CC=C", "hexa-1,3,5-triene", "己-1,3,5-三烯"),
    ("C=CCC=C", "penta-1,4-diene", "戊-1,4-二烯"),
    ("C=CC=CC", "penta-1,3-diene", "戊-1,3-二烯"),
    ("CC=CC=CC", "hexa-2,4-diene", "己-2,4-二烯"),
    ("C=CCCCC=C", "hepta-1,6-diene", "庚-1,6-二烯"),
    # negative: monoalkene / alkyne must not become polyene
    ("CC=CC", "but-2-ene", "丁-2-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", polyene__CASES)
def test_polyene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_enyne_polyyne.py
# IUPAC: P-31.1.1.1 / P-31.1.1.2
# Layer: L2,L4,L5
#
# 混合烯炔与多炔命名：ene 前 yne 后、词干插 a、多重键集合编号优先。
#
# 纯烃（alkane）、FG 段式（醇）与融合式（酸）三路均覆盖；多烯+多炔、带侧链一并。
# ==========================================================================
enyne_polyyne__CASES = [
    # 混合烯炔（纯烃）：ene 前 yne 后，双键得低位
    ("C#CCC=C", "pent-1-en-4-yne", "戊-1-烯-4-炔"),
    ("C=CC#C", "but-1-en-3-yne", "丁-1-烯-3-炔"),
    # 3 位侧基 =CH2 是外环亚甲基（P-56.4 methylidene）；旧 'methyl' 丢双键致式量错
    # （C=C(C#C)C=C 为 C6H6，饱和甲基对应 C6H10，不同分子）。
    ("C=C(C#C)C=C", "3-methylidenepent-1-en-4-yne", "3-亚甲基戊-1-烯-4-炔"),
    # 多炔：词干插 a，MULT 后缀
    ("C#CCCC#C", "hexa-1,5-diyne", "己-1,5-二炔"),
    ("C#CC#C", "buta-1,3-diyne", "丁-1,3-二炔"),
    ("C#CCC#CC", "hexa-1,4-diyne", "己-1,4-二炔"),
    # 多烯+多炔混合：dien + diyne 组合段
    ("C=C=CC#CC#C", "hepta-1,2-dien-4,6-diyne", "庚-1,2-二烯-4,6-二炔"),
    # FG 段式（醇）
    ("C#CCC(O)C=C", "hex-1-en-5-yn-3-ol", "己-1-烯-5-炔-3-醇"),
    # FG 融合式（酸）
    ("C=CC#CCC(=O)O", "hex-5-en-3-ynoic acid", "己-5-烯-3-炔酸"),
    # 多 FG（多醇/多胺）：词干走主路径，多炔/混合一并支持
    ("OCC=CCO", "but-2-ene-1,4-diol", "丁烷-2-烯-1,4-二醇"),
    ("OCC#CC#CCO", "hexa-2,4-diyne-1,6-diol", "己烷-2,4-二炔-1,6-二醇"),
    ("NCC=CCN", "but-2-ene-1,4-diamine", "丁烷-2-烯-1,4-二胺"),
    ("NCCC#CCN", "pent-2-yne-1,5-diamine", "戊烷-2-炔-1,5-二胺"),
    ("NC=CC#CCN", "pent-1-en-3-yne-1,5-diamine", "戊烷-1-烯-3-炔-1,5-二胺"),
    ("NCC=CC#CCN", "hex-2-en-4-yne-1,6-diamine", "己烷-2-烯-4-炔-1,6-二胺"),
    # 多 FG（多硫醇）：P-57 保留 e
    ("SCC=CCS", "but-2-ene-1,4-dithiol", "丁烷-2-烯-1,4-二硫醇"),
    ("SCCC#CCS", "pent-2-yne-1,5-dithiol", "戊烷-2-炔-1,5-二硫醇"),
    ("SC=CC#CCS", "pent-1-en-3-yne-1,5-dithiol", "戊烷-1-烯-3-炔-1,5-二硫醇"),
    # 对照 negative：纯烯/纯炔不得走混合路径
    ("C#C", "ethyne", "乙炔"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", enyne_polyyne__CASES)
def test_enyne_polyyne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkyne_fgs.py
# IUPAC: P-31.1 / P-63.1.1 / P-66.6.1 / P-66.5.1 / P-66.1.1 / P-14.3.4
# Layer: L4,L5
#
# Open-chain mono alkyne FG parents: alkynol / alkynal / alkynenitrile / alkynamide.
#
# Mono principal FG + mono C≡C (no C=C); FG carbon and triple-bond ends acyclic.
# Kinds stay alcohol/aldehyde/nitrile/amide; unsaturation via triple_bond + yne_locant.
#
# P-14.3.4：三碳炔的融合式后缀（-ynal/-ynenitrile/-ynamide/…）保留炔位次，
# 与对应的 prop-2-enal/prop-2-enenitrile 一致；开链烃 propyne 仍省略位次。
# ==========================================================================
alkyne_fgs__CASES = [
    # alkynols (OH principal; …-n-yn-m-ol)
    ("C#CCCO", "but-3-yn-1-ol", "丁-3-炔-1-醇"),
    ("C#CCO", "prop-2-yn-1-ol", "丙-2-炔-1-醇"),
    ("CC#CCO", "but-2-yn-1-ol", "丁-2-炔-1-醇"),
    ("C#CC(O)C", "but-3-yn-2-ol", "丁-3-炔-2-醇"),
    # alkynals
    ("C#CCC=O", "but-3-ynal", "丁-3-炔醛"),
    ("C#CC=O", "prop-2-ynal", "丙-2-炔醛"),
    # alkynenitriles
    ("C#CCC#N", "but-3-ynenitrile", "丁-3-炔腈"),
    ("C#CC#N", "prop-2-ynenitrile", "丙-2-炔腈"),
    # alkynamides
    ("C#CCC(=O)N", "but-3-ynamide", "丁-3-炔酰胺"),
    ("C#CC(=O)N", "prop-2-ynamide", "丙-2-炔酰胺"),
    # negatives: ene alcohol, sat alcohol, alkynoic acid, free alkyne
    ("C/C=C/CO", "(2E)-but-2-en-1-ol", "(2E)-丁-2-烯-1-醇"),
    ("CCCCO", "butan-1-ol", "丁-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", alkyne_fgs__CASES)
def test_alkyne_fgs(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_ring_cumulene.py
# IUPAC: P-14.4(e) / P-31.1.3 / P-22.1.2.2
# Layer: L4,L5
#
# 环内累积二烯 (ring cumulene) 编号：C=C=C 嵌入环内。
#
# 双键共享同一 sp 累积碳（如 C1=CCCCC=1 中 C0 连 (0,1) 与 (0,5) 两个 C=C），
# P-14.4(e)(ii) 编号不得让共享碳同时充当两条双键的 locant → 应为互异连续 locant
# EN cyclohexa-1,2-diene；ZH 环己-1,2-二烯。
# 回归守护：普通 1,3-/1,4- 环己二烯与单环己烯不回归。
# ==========================================================================
ring_cumulene__CASES = [
    # positive: 环内累积二烯（两个 C=C 共享 sp 中心碳，locant 须互异连续）
    ("C1=CCCCC=1", "cyclohexa-1,2-diene", "环己-1,2-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", ring_cumulene__CASES)
def test_ring_cumulene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CCCCC=1"])
def test_no_duplicate_ene_locant(smiles: str) -> None:
    """累积双键 locant 必须互异，不得退化成 1,1-diene 这类重复位次。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert "diene" in en
    assert en != "cyclohexa-1,1-diene"
    assert "1,1" not in en


# ==========================================================================
# 合并自 test_exocyclic_methylidene.py
# IUPAC: P-56.4 / P-29 (H2C= -> methylidene, preferred prefix; 不用旧名 methylene)
# Layer: L3,L5
#
# Exocyclic double bond leaves: parent=CH2 must keep the double bond (methylidene).
#
# A substituent leaf whose root atom bonds to the parent via a DOUBLE bond was
# previously collapsed to the saturated alkyl (methyl) because the anchor dummy
# bond was hardcoded single (methylcyclohexane for C=C1CCCCC1, formula loses H2).
# The exocyclic =CH2 group is IUPAC 'methylidene' (P-56.4); endocyclic/chain
# alkenes and saturated twins must stay unchanged.
# ==========================================================================
exocyclic_methylidene__POSITIVE = [
    # exocyclic =CH2 on a ring -> methylidene retained
    ("C=C1CCCCC1", "methylidenecyclohexane", "亚甲基环己烷"),
    # exocyclic =CH2 on a chain near FG + stereo
    ("C=C(C[C@H](N)C(=O)O)C(=O)O",
     "(2S)-2-amino-4-methylidenepentanedioic acid",
     "(2S)-2-氨基-4-亚甲基戊二酸"),
    # nested: =CH2 on a cyclopropyl substituent (recursive path)
    ("C=C1CC1CC(=O)C(=O)O",
     "3-(2-methylidenecyclopropyl)-2-oxopropanoic acid",
     "3-(2-亚甲基环丙基)-2-氧代丙酸"),
]

exocyclic_methylidene__NEGATIVE = [
    # saturated twin must NOT take the ylidene route
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    # endocyclic/chain alkenes stay alkene names
    ("C1=CCCCC1", "cyclohexene", "环己烯"),
    ("C1=CCCC1", "cyclopentene", "环戊烯"),
    ("C=C(C)C", "2-methylprop-1-ene", "2-甲基丙-1-烯"),
]


@pytest.mark.parametrize("smiles,en,zh", exocyclic_methylidene__POSITIVE)
def test_exocyclic_methylidene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", exocyclic_methylidene__NEGATIVE)
def test_neighbors_not_ylidene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_exocyclic_not_saturated() -> None:
    """Guard: exocyclic =CH2 must not be reported as the saturated methyl twin."""
    r = SMILESNNamer().name("C=C1CCCCC1")
    assert normalize_en(r.en) != normalize_en("methylcyclohexane")
    assert "methylidene" in normalize_en(r.en)


# ==========================================================================
# 合并自 test_long_chain_stems.py
# IUPAC: P-14.2.1
# Layer: L5
#
# Long-chain open parent stems C11–C35 (numerical multipliers + alkane stems).
#
# Extends alkane / alcohol / acid / aldehyde / haloalkane naming past C10 so
# L5 no longer returns empty for unsupported n_carbons>10. Chinese multi-char
# stems (十一…三十五) must not be truncated via zh[0].
# ==========================================================================
long_chain_stems__CASES = [
    ("CCCCCCCCCCCCCCI", "1-iodotetradecane", "1-碘十四烷"),
    ("C=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "pentatriacont-1-ene", "三十五-1-烯"),
    ("CCCCCCCCCCC", "undecane", "十一烷"),
    ("CCCCCCCCCCC(=O)O", "undecanoic acid", "十一酸"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=O", "triacontanal", "三十醛"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", "tetratriacontan-1-ol", "三十四-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", long_chain_stems__CASES)
def test_long_chain_stems(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
