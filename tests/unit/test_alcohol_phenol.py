# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_alkanediol.py: Unsubstituted open-chain saturated alkanediols (exactly two OH).
test_alkanetriol.py: Unsubstituted open-chain saturated alkanetriols (exactly three OH).
test_alkenol.py: Open-chain unsaturated monoalcohols (alkenols / polyalkenols).
test_alkanethiol.py: Unsubstituted open-chain monovalent thiols (alkanethiols).
test_benzenediol.py: Unsubstituted benzene-a,b-diol (exactly two phenolic OH on benzene).
test_phenol_aniline.py: Retained parent names phenol / aniline (aromatic mono-OH / mono-NH2).
test_aminophenol.py: Aminophenols: phenol parent + amino prefix (OH senior to NH2).
test_ring_ether_alpha_ol.py: 环醚 α-羟基(半缩醛型)醇:醇碳同时邻接另一单键氧(环醚/醚 O)时,醇 OH 氧不得游离出
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_alkanediol.py
# IUPAC: P-63.1.1
# Layer: L2,L4,L5
#
# Unsubstituted open-chain saturated alkanediols (exactly two OH).
#
# Parent chain through both OH carbons; name alkane-{a},{b}-diol /
# {烷}-{a},{b}-二醇. Side chains use existing L3 alkyl prefixes.
# ==========================================================================
alkanediol__CASES = [
    ("CC(O)C(C)O", "butane-2,3-diol", "丁烷-2,3-二醇"),
    ("OCC(C)O", "propane-1,2-diol", "丙烷-1,2-二醇"),
    ("CC(C)C(O)CO", "3-methylbutane-1,2-diol", "3-甲基丁烷-1,2-二醇"),
    ("OCCCCO", "butane-1,4-diol", "丁烷-1,4-二醇"),
    ("CC(O)C", "propan-2-ol", "丙-2-醇"),
    ("CC(C)(O)O", "propane-2,2-diol", "丙烷-2,2-二醇"),  # 同碳二醇：位次须重复写出
    ("OC1(O)CCCCC1", "cyclohexane-1,1-diol", "环己烷-1,1-二醇"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanediol__CASES)
def test_alkanediol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanetriol.py
# IUPAC: P-63.1.1
# Layer: L2,L4,L5
#
# Unsubstituted open-chain saturated alkanetriols (exactly three OH).
#
# Parent chain through all three OH carbons; name alkane-{a},{b},{c}-triol /
# {烷}-{a},{b},{c}-三醇. Diols and monoalcohols must not break.
# ==========================================================================
alkanetriol__CASES = [
    ("OCC(O)CCCO", "pentane-1,2,5-triol", "戊烷-1,2,5-三醇"),
    ("CC(O)C(O)CO", "butane-1,2,3-triol", "丁烷-1,2,3-三醇"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanetriol__CASES)
def test_alkanetriol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenol.py
# IUPAC: P-63.1.1 / P-31.1 / P-93.4
# Layer: L2,L4,L5
#
# Open-chain unsaturated monoalcohols (alkenols / polyalkenols).
#
# Alcohol is the principal characteristic group; non-aromatic C=C is inserted as
# -ene_locant-en-OH_locant-ol (mono) or -a,b-dien-OH-ol (poly). Parent chain covers
# the OH carbon and all double-bond carbons. Numbering: lowest OH locant first,
# then lowest ene set. BondStereo → (E)-/(Z)- or (2E,6Z)- prefixes (P-93.4).
# ==========================================================================
alkenol__CASES = [
    # C2 单烯醇：烯 1-2 + OH 1 无歧义，位次省略融合 ethenol (P-14.3.4)
    ("C=CO", "ethenol", "乙烯醇"),
    ("C=CCO", "prop-2-en-1-ol", "丙-2-烯-1-醇"),
    ("CC(O)C=C", "but-3-en-2-ol", "丁-3-烯-2-醇"),
    ("C=CCCCCO", "hex-5-en-1-ol", "己-5-烯-1-醇"),
    ("CC(C)=CCO", "3-methylbut-2-en-1-ol", "3-甲基丁-2-烯-1-醇"),
    # E/Z mono
    (r"C(\C=C\CCCCCCCCC)O", "(2E)-dodec-2-en-1-ol", "(2E)-十二-2-烯-1-醇"),
    (r"CCCC/C=C\CCCCCCCCCCCCO", "(13Z)-octadec-13-en-1-ol", None),
    (r"C/C=C/CCO", "(3E)-pent-3-en-1-ol", "(3E)-戊-3-烯-1-醇"),
    (r"CC/C=C\CC/C=C/CO", "(2E,6Z)-nona-2,6-dien-1-ol", "(2E,6Z)-壬-2,6-二烯-1-醇"),
    (
        r"CC/C=C\C/C=C\C/C=C\CCCCCCCCO",
        "(9Z,12Z,15Z)-octadeca-9,12,15-trien-1-ol",
        None,
    ),
    # polyalkenols beyond triene (multiplicity-generic mult-seg)
    ("C=CC=CC=CC=CCO", "nona-2,4,6,8-tetraen-1-ol", "壬-2,4,6,8-四烯-1-醇"),
    ("C=CC=CC=CC=CC=CCO", "undeca-2,4,6,8,10-pentaen-1-ol", "十一-2,4,6,8,10-五烯-1-醇"),
    ("C=CC=CC=CC", "hepta-1,3,5-triene", "庚-1,3,5-三烯"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenol__CASES)
def test_alkenol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanethiol.py
# IUPAC: P-63.1.5
# Layer: L1,L2,L4,L5
#
# Unsubstituted open-chain monovalent thiols (alkanethiols).
#
# Exactly one –SH (S with one C neighbor and ≥1 H); no other principal FG;
# saturated acyclic. Parent chain through SH-attached carbon; locants like alcohol.
# C1–C2 omit locant (methanethiol/ethanethiol); C≥3: alkane-n-thiol / 首字-n-硫醇.
# ==========================================================================
alkanethiol__CASES = [
    # C2 烯硫醇：烯 1-2 + SH 1 无歧义，位次省略融合 ethenethiol (P-14.3.4)
    ("C=CS", "ethenethiol", "乙烯硫醇"),
    # alkanethiol 正例（P-63.1.5）：C1–C2 omit locant，C≥3 用 alkane-n-thiol
    ("CS", "methanethiol", "甲硫醇"),
    ("CCS", "ethanethiol", "乙硫醇"),
    ("CCCS", "propane-1-thiol", "丙-1-硫醇"),
    ("CCCCCS", "pentane-1-thiol", "戊-1-硫醇"),
    ("CC(C)S", "propane-2-thiol", "丙-2-硫醇"),
    ("CCCC(CC)S", "hexane-3-thiol", "己-3-硫醇"),
    # 二硫醇（数量后缀生成式，对齐 amine）
    ("SCCS", "ethane-1,2-dithiol", "乙烷-1,2-二硫醇"),
    ("SCCCCS", "butane-1,4-dithiol", "丁烷-1,4-二硫醇"),
    ("CC(C)(S)S", "propane-2,2-dithiol", "丙烷-2,2-二硫醇"),  # 同碳二硫醇：位次须重复写出
]


@pytest.mark.parametrize("smiles,en,zh", alkanethiol__CASES)
def test_alkanethiol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzenediol.py
# IUPAC: P-63.1.2
# Layer: L2,L3,L4,L5
#
# Unsubstituted benzene-a,b-diol (exactly two phenolic OH on benzene).
#
# Systematic names only: benzene-1,2-diol / benzene-1,3-diol / benzene-1,4-diol.
# Chinese uses 二酚 (gold), not 二醇. No retained catechol/resorcinol/hydroquinone.
# ==========================================================================
benzenediol__CASES = [
    # positive: benzene-diol/triol 走主路径 _RING_STEM 词干 + aromatic 酚（非 variant）
    ("Oc1ccccc1O", "benzene-1,2-diol", "苯-1,2-二酚"),
    ("Oc1cc(O)cc(O)c1", "benzene-1,3,5-triol", "苯-1,3,5-三酚"),
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
]


@pytest.mark.parametrize("smiles,en,zh", benzenediol__CASES)
def test_benzenediol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_phenol_aniline.py
# IUPAC: P-63.1.4 / P-62.2.1.1.1
# Layer: L2,L4,L5
#
# Retained parent names phenol / aniline (aromatic mono-OH / mono-NH2).
#
# Parent = benzene with principal characteristic OH (phenol) or NH2 (aniline).
# Allow 0–2 extra ring halo and/or methyl substituents; FG fixed as locant 1.
# ==========================================================================
phenol_aniline__CASES = [


    # positive: monohalo / monomethyl
    ("Oc1ccc(Cl)cc1", "4-chlorophenol", "4-氯苯酚"),
    ("Cc1ccccc1O", "2-methylphenol", "2-甲基苯酚"),
    ("Nc1ccc(C)cc1", "4-methylaniline", "4-甲基苯胺"),
    # positive: dimethyl phenol
    ("Cc1cc(C)cc(O)c1", "3,5-dimethylphenol", "3,5-二甲基苯酚"),
]


@pytest.mark.parametrize("smiles,en,zh", phenol_aniline__CASES)
def test_phenol_aniline(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_phenol_not_chain_alcohol() -> None:
    """Aromatic OH must not expand ring into chain alcohol parent."""
    r = SMILESNNamer().name("Oc1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "phenol"
    assert "hexanol" not in en
    assert "nonan" not in en


# ==========================================================================
# 合并自 test_aminophenol.py
# IUPAC: P-63.1.4 / P-62.5
# Layer: L2,L3,L4,L5
#
# Aminophenols: phenol parent + amino prefix (OH senior to NH2).
#
# Benzene with exactly one phenolic OH and one primary ring NH2.
# OH is locant 1; amino gets 2/3/4. Optional extra halo/methyl/nitro
# still limited by arene FG sub cap (amino counts as one extra).
# ==========================================================================
aminophenol__CASES = [
    # positive: unsubstituted aminophenols (o/m/p)
    ("Nc1ccc(O)cc1", "4-aminophenol", "4-氨基苯酚"),
    ("Nc1ccccc1O", "2-aminophenol", None),
    ("Nc1cc(O)ccc1", "3-aminophenol", None),





]


@pytest.mark.parametrize("smiles,en,zh", aminophenol__CASES)
def test_aminophenol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_ring_ether_alpha_ol.py
# IUPAC: P-63.1.1 / P-22.1.1
# Layer: L2 (parent_ownership) + L5
#
# 环醚 α-羟基(半缩醛型)醇:醇碳同时邻接另一单键氧(环醚/醚 O)时,醇 OH 氧不得游离出
# owned_atoms 被 coverage 补齐阶段命名成假 hydroxy 前缀(与 -ol 后缀双算)。根因:_single_o_idx
# 盲取醇碳首个单键 O,歧义时取到醚氧。
# ==========================================================================
ring_ether_alpha_ol__CASES = [
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


@pytest.mark.parametrize("smiles,en,zh", ring_ether_alpha_ol__CASES)
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
        ("OCC(O)C", "propane-1,2-diol", "丙烷-1,2-二醇"),
    ],
)
def test_unaffected_alcohols_still_fine(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
