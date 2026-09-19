# 合并自 10 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_amine.py: Simple acyclic primary monoamines (alkanamines).
test_polyamine.py: Unsubstituted open-chain saturated polyamines (triamine, tetraamine).
test_alkanediamine.py: Unsubstituted open-chain saturated alkanediamines (exactly two primary amines).
test_alkenamine.py: Unsaturated primary amines (alkenamines).
test_sec_amine.py: Open-chain secondary monoamine scope (P-62.2.2.1) — negative guard only.
test_tert_amine.py: Open-chain tertiary monoamine scope (P-62.2.2.1) — negative guard only.
test_tetraalkylammonium.py:
test_zh_amino_bridge_root.py: 桥后缀（氨基/氧基/硫基）前中文烃基名的「基」字去留。
test_arylalkyl_amine.py: Sec-amine: N-端取代基命名（简单 2° 胺）与芳烷基对照。
test_benzene_poly_diamine.py: Multi-substituted benzene (cap 4) and benzene-1,n-diamine parents.
"""
from __future__ import annotations

import pytest

from namepredict.layer1.analyzer import analyze
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG, inventory_from_info
from namepredict.layer2.parent_select import _collect_candidates
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# 合并自 test_mono_amine.py
# IUPAC: P-62.2.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic primary monoamines (alkanamines).
#
# Primary amine: N with exactly one carbon neighbor and ≥2 H; not amide nitrogen.
# Parent chain through the carbon attached to N; kind=amine.
# C1/C2 omit locant (methanamine/ethanamine); C≥3 use …-n-amine / …-n-胺.
# ==========================================================================
mono_amine__CASES = [
    # positive: primary monoamines (straight + optional branched)
    ("CN", "methanamine", "甲胺"),
    ("CCCCCN", "pentan-1-amine", "戊-1-胺"),
    ("CC(C)N", "propan-2-amine", "丙-2-胺"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_amine__CASES)
def test_mono_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_amide_not_named_as_amine() -> None:
    """Amide nitrogen (N next to carbonyl) must not become *amine."""
    r = SMILESNNamer().name("CC(=O)N")
    assert r.success
    assert "amine" not in normalize_en(r.en)


# ==========================================================================
# 合并自 test_polyamine.py
# IUPAC: P-62.2.1
# Layer: L2,L3,L4,L5
#
# Unsubstituted open-chain saturated polyamines (triamine, tetraamine).
#
# Parent chain through all amine carbons; name alkane-{a},{b},{c}-triamine /
# {烷}-{a},{b},{c}-三胺 (triamine) or alkane-{a},{b},{c},{d}-tetraamine /
# {烷}-{a},{b},{c},{d}-四胺 (tetraamine).
#
# Diamines and monoamines must not break.
# ==========================================================================
polyamine__CASES = [
    # triamine: unsubstituted open-chain
    ("NCC(N)CN", "propane-1,2,3-triamine", "丙烷-1,2,3-三胺"),
    ("NCCC(N)CCN", "pentane-1,3,5-triamine", "戊烷-1,3,5-三胺"),
    ("NCC(N)(C)CN", "2-methylpropane-1,2,3-triamine", "2-甲基丙烷-1,2,3-三胺"),
    ("NCC(N)C(N)C", "butane-1,2,3-triamine", "丁烷-1,2,3-三胺"),
    # multi-substituted triamine edge cases
    ("NCC(N)(CC)CN", "2-ethylpropane-1,2,3-triamine", "2-乙基丙烷-1,2,3-三胺"),
    ("NCC(C)C(N)CN", "3-methylbutane-1,2,4-triamine", "3-甲基丁烷-1,2,4-三胺"),
]


@pytest.mark.parametrize("smiles,en,zh", polyamine__CASES)
def test_polyamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanediamine.py
# IUPAC: P-62.2.1
# Layer: L2,L4,L5
#
# Unsubstituted open-chain saturated alkanediamines (exactly two primary amines).
#
# Parent chain through both amine carbons; name alkane-{a},{b}-diamine /
# {烷}-{a},{b}-二胺. Side chains use existing L3 alkyl prefixes.
# ==========================================================================
alkanediamine__CASES = [
    ("NCCCN", "propane-1,3-diamine", "丙烷-1,3-二胺"),
    ("NCCCCCN", "pentane-1,5-diamine", "戊烷-1,5-二胺"),
    ("NC(C)CCN", "butane-1,3-diamine", "丁烷-1,3-二胺"),
    ("NCC(C)CN", "2-methylpropane-1,3-diamine", "2-甲基丙烷-1,3-二胺"),
    ("CC(N)C", "propan-2-amine", "丙-2-胺"),
    ("C1CCC(N)CC1", "cyclohexanamine", "环己胺"),
    ("CC(C)(N)N", "propane-2,2-diamine", "丙烷-2,2-二胺"),  # 同碳二胺：位次须重复写出
]


@pytest.mark.parametrize("smiles,en,zh", alkanediamine__CASES)
def test_alkanediamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenamine.py
# IUPAC: P-62.2 / P-31.1
# Layer: L2,L4,L5
#
# Unsaturated primary amines (alkenamines).
#
# Amines use the segment-style ene insertion (same as alcohols):
# but-3-en-1-amine. Previously the amine entry had no ene segment, so double
# bonds were silently dropped (butan-1-amine).
# ==========================================================================
alkenamine__CASES = [
    # C≥3：烯/胺位次需区分位置异构 → 保留 (P-62.2.6.2)
    ("C=CCCN", "but-3-en-1-amine", "丁-3-烯-1-胺"),

    # C2：烯只能 1-2、胺后缀锚定 1，位次无歧义省略融合 (P-14.3.4)
    ("C=CN", "ethenamine", "乙烯胺"),
    # #502：短链省略仅去掉两个 '1'，另一端硝基位次 2 仍保留
    (r"C(C1=CC=CC=C1)N/C=C\[N+](=O)[O-]",
     "(1Z)-N-benzyl-2-nitroethenamine", "(1Z)-N-苄基-2-硝基乙烯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenamine__CASES)
def test_alkenamine(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_sec_amine.py
# IUPAC: P-62.2.2.1
# Layer: L1,L2,L3,L5
#
# Open-chain secondary monoamine scope (P-62.2.2.1) — negative guard only.
#
# Positive secondary-amine cases are not covered in this file. The retained
# cases assert that a primary monoamine (CCN) and a diamine (NCCN) must not be
# named as N-alkylalkanamines.
# ==========================================================================
sec_amine__CASES = [
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", sec_amine__CASES)
def test_sec_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_tert_amine.py
# IUPAC: P-62.2.2.1
# Layer: L1,L2,L3,L5
#
# Open-chain tertiary monoamine scope (P-62.2.2.1) — negative guard only.
#
# Positive tertiary-amine cases are not covered in this file. The retained
# cases assert that a primary monoamine (CCN) and a diamine (NCCN) must not be
# named as N,N-dialkylalkanamines.
# ==========================================================================
tert_amine__CASES = [
    # negative: primary / secondary / diamine must not break
    ("CCN", "ethanamine", "乙胺"),
]


@pytest.mark.parametrize("smiles,en,zh", tert_amine__CASES)
def test_tert_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_tetraalkylammonium.py
# IUPAC: P-41 cation / tetraalkylammonium
# Layer: L1,L2,L5
# ==========================================================================
def test_quaternary_is_not_neutral_amine():
    info = analyze(Chem.MolFromSmiles("C[N+](C)(C)C"))
    assert not inventory_from_info(info).occurrences(FG.AMINE)


# ==========================================================================
# 合并自 test_zh_amino_bridge_root.py
# IUPAC: P-62.2 / P-15.4.1
# Layer: L5
#
# 桥后缀（氨基/氧基/硫基）前中文烃基名的「基」字去留。
#
# 规则：桥后缀直接相连的烃基名省略尾「基」字（甲基→甲、环己基→环己、
# 丙-2-基→丙-2-），与 tiers gold 口径一致。母体为胺时 N-取代前缀
# （N,N-二甲基苯胺）不受此规则影响，仍保留「基」。
# ==========================================================================
zh_amino_bridge_root__CASES = [
    # ── 正例：N-多取代简单烷氨基作为取代基，中文省「基」 ──
    ("CN(C)CCO", "2-(dimethylamino)ethanol", "2-(二甲氨基)乙醇"),
    ("CCCN(CCC)c1ccccc1C(=O)O", "2-(dipropylamino)benzoic acid", "2-(二丙氨基)苯甲酸"),
    # ── 近邻：位次化/环系烃基名同样省「基」（tiers gold：丙-2-氧基、环己氧基、环丙硫基） ──
    ("CC(C)NCCO", "2-(propan-2-ylamino)ethanol", "2-(丙-2-氨基)乙醇"),
    ("C1CCCCC1NCCO", "2-(cyclohexylamino)ethanol", "2-(环己氨基)乙醇"),
    # ── 近邻负例：胺为母体时的 N-取代前缀仍保留「基」 ──
    ("CN(C)c1ccccc1", "N,N-dimethylaniline", "N,N-二甲基苯胺"),
    ("CN(C)CCN(C)C", "N,N,N',N'-tetramethylethane-1,2-diamine",
     "N,N,N',N'-四甲基乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", zh_amino_bridge_root__CASES)
def test_zh_bridge_root(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_arylalkyl_amine.py
# IUPAC: P-62.2.2.1 / P-29.3
# Layer: L2,L3
#
# Sec-amine: N-端取代基命名（简单 2° 胺）与芳烷基对照。
#
# 含芳基臂的二级胺（N-ethyl-1-phenylethanamine 等）通过 L2 P-45.2.1 流水线
# 选出前缀取代基团数目更多（芳环臂作为母体骨架取代基）的母体链，主链臂
# 优先保留为母体、另一臂作 N-取代基。
# ==========================================================================
arylalkyl_amine__CASES = [
    # nitriles from user set (should already work)
    (
        "C(C)C=1C=C(C=CC1)CC#N",
        "2-(3-ethylphenyl)acetonitrile",
        "2-(3-乙基苯基)乙腈",
    ),
    (
        "ClC1=C(C(=CC=C1)F)CC#N",
        "2-(2-chloro-6-fluorophenyl)acetonitrile",
        "2-(2-氯-6-氟苯基)乙腈",
    ),
    # ester: prefix inserts between alkyl and acyl (user ex7)
    (
        "COC(C(CCC1=CC=CC=C1)=O)=O",
        "methyl 2-oxo-4-phenylbutanoate",
        "2-氧代-4-苯基丁酸甲酯",
    ),
    # 二级胺：N- 取代基前缀
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("NCCc1ccccc1", "2-phenylethanamine", "2-苯基乙胺"),
    # 含芳基臂的二级胺：母体选前缀取代基更多的链（P-45.2.1），芳环臂作母体骨架取代基
    (
        "C(C)NC(C)C1=CC=C(C=C1)OC",
        "N-ethyl-1-(4-methoxyphenyl)ethanamine",
        "N-乙基-1-(4-甲氧基苯基)乙胺",
    ),
    (
        "C(C)NC(C)c1ccccc1",
        "N-ethyl-1-phenylethanamine",
        "N-乙基-1-苯基乙胺",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", arylalkyl_amine__CASES)
def test_arylalkyl_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzene_poly_diamine.py
# IUPAC: P-14.3.4 / P-62.2.1 / P-61.5
# Layer: L2,L4,L5
#
# Multi-substituted benzene (cap 4) and benzene-1,n-diamine parents.
#
# Simple ring subs (halo / alkyl / alkoxy / nitro / CF3) up to 4 on benzene.
# Two primary amines on benzene → benzene-a,b-diamine (≤2 simple halo/methyl).
# Gold chinese_name wins when it conflicts with pure system zh (anisole-style).
# ==========================================================================
benzene_poly_diamine__CASES = [
    # positive: 4-sub benzene
    (
        "ClC1=C(C(=CC(=C1)OC)[N+](=O)[O-])C",
        "1-chloro-5-methoxy-2-methyl-3-nitrobenzene",
        "1-氯-5-甲氧基-2-甲基-3-硝基苯",
    ),
    (
        "FC1=C(C=C(C(=C1)F)[N+](=O)[O-])[N+](=O)[O-]",
        "1,5-difluoro-2,4-dinitrobenzene",
        "1,5-二氟-2,4-二硝基苯",
    ),
    # negative: keep existing correct behaviour
    ("Clc1ccc(Cl)c(Cl)c1", "1,2,4-trichlorobenzene", "1,2,4-三氯苯"),
]


@pytest.mark.parametrize("smiles,en,zh", benzene_poly_diamine__CASES)
def test_benzene_poly_diamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
