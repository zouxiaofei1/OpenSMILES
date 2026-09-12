# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_nitrile.py: Simple acyclic mononitriles (alkanenitriles).
test_alkenenitrile.py: Open-chain monounsaturated mononitriles (alkenenitriles).
test_nitrile_cyano_demote.py: 非主官能团的腈须前缀化为 cyano，不得被词干链吞碳把 N 悬空误命名成 amino。
test_cyano_as_leaf.py: 端位 C≡N 在更高优先级主基团（自由基/羧酸/酰胺/酯…）下退为 cyano 前缀叶：
test_isocyanate.py: Isocyanate / isothiocyanate (R–N=C=X, X=O/S): functional class + arene prefix.
test_guanidine.py: Guanidine functional parent H2N–C(=NH)–NH2 (retained name guanidine).
test_hydrazine.py: Hydrazine scope (P-68.3.1.2) — negative guard only.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mono_nitrile.py
# IUPAC: P-66.5.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic mononitriles (alkanenitriles).
#
# Nitrile carbon is C≡N; chain through that C. Retained acetonitrile;
# C≥3 systematic …nitrile / …腈. Not alkyne (C≡C).
# ==========================================================================
mono_nitrile__CASES = [
    # positive
    ("C#N", "formonitrile", "甲腈"),
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("CCCCC#N", "pentanenitrile", "戊腈"),
    ("CCCCCC#N", "hexanenitrile", "己腈"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_nitrile__CASES)
def test_mono_nitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenenitrile.py
# IUPAC: P-66.5.1 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated mononitriles (alkenenitriles).
#
# Nitrile is the principal characteristic group (C≡N carbon = locant 1); one
# non-aromatic C=C is expressed as -n-enenitrile / -n-烯腈 with the lower
# double-bond carbon locant. No (E)/(Z) stereodescriptors this cycle.
# ==========================================================================
alkenenitrile__CASES = [
    # positive: acyclic mono-alkenenitriles (no E/Z)
    ("C=CC#N", "prop-2-enenitrile", "丙-2-烯腈"),
    ("CC=CC#N", "but-2-enenitrile", "丁-2-烯腈"),
    ("C=CCC#N", "but-3-enenitrile", "丁-3-烯腈"),
    ("CCC=CC#N", "pent-2-enenitrile", "戊-2-烯腈"),
    ("C=CCCC#N", "pent-4-enenitrile", "戊-4-烯腈"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenenitrile__CASES)
def test_alkenenitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_nitrile_cyano_demote.py
# IUPAC: P-61.1.3,P-66.1.1
# Layer: L1,L2,L3
#
# 非主官能团的腈须前缀化为 cyano，不得被词干链吞碳把 N 悬空误命名成 amino。
#
# 腈碳(C≡N)是「自带碳的前缀叶」：当自由基(p41=1)/羧酸/酰胺等更高优先级主基团存在时，
# 腈退出主基团并保留 cyano 叶身份——其碳须从开链主链剔除（同中性 COOH 的 carboxy 叶，
# P-61.1.3），否则如 `N#CCC*`(→*CH₂CH₂CN) 会被读成饱和丙基链 + 端 N 而错名 3-aminopropyl；
# `HOOC-CH₂CN` 也被错名成 3-aminopropanoic acid。
#
# - 链端/支链腈同规则：`*CH₂CH₂CN` → 2-cyanoethyl；`HOOC-CH₂CN` → 2-cyanoacetic acid；
#   `HOOC-CH(CN)CH₃` → 2-cyanopropanoic acid（腈碳不得占主链位）。
# - 无更高优先级基团时腈仍是主基团后缀（butanenitrile/…），不得 cyano 化。
# ==========================================================================
nitrile_cyano_demote__POS_FRAG = [
    ("N#CCC*", "2-cyanoethyl", "2-氰基乙基"),
    ("N#CC*", "cyanomethyl", "氰基甲基"),
    ("N#CCCC*", "3-cyanopropyl", "3-氰基丙基"),
]


@pytest.mark.parametrize("smiles,en,zh", nitrile_cyano_demote__POS_FRAG)
def test_cyano_radical_fragment(smiles: str, en: str, zh: str) -> None:
    """锚定自由基残基：腈作前缀叶，N 不得误命名成 amino。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh


# 整分子正例：酸/酰胺主基团压制腈（链端或支链），腈碳不进主链。
nitrile_cyano_demote__POS_WHOLE = [
    ("OC(=O)CC#N", "2-cyanoacetic acid", "2-氰基乙酸"),
    ("OC(=O)CCC#N", "3-cyanopropanoic acid", "3-氰基丙酸"),
    ("CC(C#N)C(=O)O", "2-cyanopropanoic acid", "2-氰基丙酸"),
    ("OC(=O)C(C#N)CC", "2-cyanobutanoic acid", "2-氰基丁酸"),
    ("NC(=O)CC#N", "2-cyanoacetamide", "2-氰基乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", nitrile_cyano_demote__POS_WHOLE)
def test_cyano_whole_molecule(smiles: str, en: str, zh: str) -> None:
    """整分子：腈被更高优先级基团压制 → cyano 前缀叶，主基团后缀不变。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh


# 负例：无更高优先级基团时腈保持主基团（后缀 nitrile），不得 cyano 化 / amino 化。
nitrile_cyano_demote__NEG = [
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("NCC#N", "aminoacetonitrile", "氨基乙腈"),
    ("O=CCC#N", "3-oxopropanenitrile", "3-氧代丙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", nitrile_cyano_demote__NEG)
def test_nitrile_principal_untouched(smiles: str, en: str, zh: str) -> None:
    """腈作主基团（后缀 -nitrile/-腈）时不受 cyano 前缀化影响。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert r.zh == zh


# ==========================================================================
# 合并自 test_cyano_as_leaf.py
# IUPAC: P-41（HOOC-CH₂-CH₂- → 2-carboxyethyl 类比）/ P-61.1.3
# Layer: L1,L2,L3
#
# 端位 C≡N 在更高优先级主基团（自由基/羧酸/酰胺/酯…）下退为 cyano 前缀叶：
#
# 词干链不得吞掉腈碳——否则开链骨架把腈碳当饱和烷基碳、N 悬空误命名成 amino
# （N#CCC* 曾错成 3-aminopropyl，应为 2-cyanoethyl）。修复点：L1
# _arbitrate_parts 腈降级进 demoted_nitriles（同羧酸叶），L2 parent_skeleton
# 把降级腈碳剔除出开链主链。腈自身为主基团（无更高类）时不受影响。
# ==========================================================================
cyano_as_leaf__POS_FRAG = [
    ("N#CC*", "cyanomethyl", "氰基甲基"),
    ("N#CCC*", "2-cyanoethyl", "2-氰基乙基"),
    ("N#CCCC*", "3-cyanopropyl", "3-氰基丙基"),
    ("N#CCCCC*", "4-cyanobutyl", "4-氰基丁基"),
]


@pytest.mark.parametrize("smiles,en,zh", cyano_as_leaf__POS_FRAG)
def test_cyano_alkyl_fragment(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子：腈作 cyano 前缀挂在酸/酰胺/酯主基上。
cyano_as_leaf__POS_WHOLE = [
    ("N#CCC(=O)O", "2-cyanoacetic acid", "2-氰基乙酸"),
    ("N#CCCC(=O)O", "3-cyanopropanoic acid", "3-氰基丙酸"),
    ("NC(=O)CC#N", "2-cyanoacetamide", "2-氰基乙酰胺"),
    ("N#CCC(=O)OC", "methyl cyanoacetate", "氰基乙酸甲酯"),
    ("N#CCCc1ccccc1C(=O)O", "2-(2-cyanoethyl)benzoic acid", "2-(2-氰基乙基)苯甲酸"),
    ("N#CCCc1ccc(C(=O)O)cc1", "4-(2-cyanoethyl)benzoic acid", "4-(2-氰基乙基)苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", cyano_as_leaf__POS_WHOLE)
def test_cyano_leaf_whole_molecule(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 近邻负例：腈仍是主基团（无更高类）或芳环 cyano 前缀，不得误伤。
cyano_as_leaf__NEG = [
    # 醇比腈低优先，腈保持主基。
    ("N#CCO", "hydroxyacetonitrile", "羟基乙腈"),
    # 醛（p41 高于腈）压成 oxo，腈保持主基。
    ("N#CCCCC=O", "5-oxopentanenitrile", "5-氧代戊腈"),
    # 芳环 cyano 前缀（酸主基）不变。
    ("N#Cc1ccc(C(=O)O)cc1", "4-cyanobenzoic acid", "4-氰基苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", cyano_as_leaf__NEG)
def test_cyano_not_demoted(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_isocyanate.py
# IUPAC: P-61.9
# Layer: L1,L2,L3,L5
#
# Isocyanate / isothiocyanate (R–N=C=X, X=O/S): functional class + arene prefix.
#
# Alkyl mono: functional-class '{alkyl} isocyanate/isothiocyanate'.
# Aryl: isocyanato/isothiocyanato prefixes on benzene (with halo/methyl).
# Must not mis-detect as amide (N=C=O) or collapse isothiocyanate to alkane.
# ==========================================================================
isocyanate__CASES = [
    # positive: alkyl isocyanates (functional class)
    ("CCCCN=C=O", "butyl isocyanate", "异氰酸丁酯"),
    ("CN=C=O", "methyl isocyanate", "异氰酸甲酯"),
    # positive: aryl isocyanate / isothiocyanate
    ("c1ccccc1N=C=O", "isocyanatobenzene", "异氰酸苯酯"),
    ("c1ccccc1N=C=S", "isothiocyanatobenzene", "异硫氰酸苯酯"),
    # positive: fused iso keeps relative locants (not mono-benzene omit)
    (
        "FC(C1=CC=C(C=C1)N=C=S)(F)F",
        "4-(trifluoromethyl)isothiocyanatobenzene",
        "4-三氟甲基异硫氰酸苯酯",
    ),
    ("Cc1ccc(N=C=S)cc1", "4-methylisothiocyanatobenzene", "4-甲基异硫氰酸苯酯"),
    # positive: multi-sub arene (prefix style; locants match gold EN)
    (
        "CC1=C(C=CC(=C1)C)N=C=S",
        "2,4-dimethylisothiocyanatobenzene",
        "2,4-二甲基异硫氰酸苯酯",
    ),
    (
        "ClC1=C(C=C(C=C1)C)N=C=O",
        "1-chloro-2-isocyanato-4-methylbenzene",
        "1-氯-2-异氰酸根合-4-甲基苯",
    ),
    # negative: true amide / nitrile / nitro must not regress
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCCNC=O", "N-propylformamide", "N-丙基甲酰胺"),
    ("c1ccccc1C#N", "benzonitrile", "苯甲腈"),
    ("c1ccccc1[N+](=O)[O-]", "nitrobenzene", "硝基苯"),
]


# ==========================================================================
# 合并自 test_guanidine.py
# IUPAC: P-66.4.1.2.1
# Layer: L1,L2,L5
#
# Guanidine functional parent H2N–C(=NH)–NH2 (retained name guanidine).
#
# Scope (first cut):
# - unsubstituted guanidine
# - mono N-aryl guanidine (unfused Ph, Me/halo ≤2) → 1-phenylguanidine style
# - N-arylsulfonyl guanidine Ar–SO2–NH–C(=NH)NH2 → 1-(…sulfonyl)guanidine dual
# Must not regress methanamine, aniline, urea, benzenesulfonamide, sulfonyl chloride.
# ==========================================================================
guanidine__CASES = [
    # positive: unsubstituted
    ("N=C(N)N", "guanidine", "胍"),
    # positive: mono N-aryl
    ("N=C(N)Nc1ccccc1", "1-phenylguanidine", "1-苯基胍"),
    # positive: N-arylsulfonyl (benchmark dual)
    (
        "ClC1=C(C=CC=C1)S(=O)(=O)NC(=N)N",
        "1-(2-chlorophenylsulfonyl)guanidine",
        "1-(2-氯苯基磺酰基)胍",
    ),
    # negative: primary amine / aniline / urea / sulfonamide / sulfonyl chloride
    ("CN", "methanamine", "甲胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
    ("NC(=O)N", "urea", "脲"),
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride", "甲磺酰氯"),
]


# ==========================================================================
# 合并自 test_hydrazine.py
# IUPAC: P-68.3.1.2
# Layer: L0,L1,L2,L5
#
# Hydrazine scope (P-68.3.1.2) — negative guard only.
#
# Positive hydrazine cases are not covered in this file. The retained cases
# assert that acetamide and aniline are not named as hydrazines.
# ==========================================================================
hydrazine__CASES = [
    # negative: amide / aniline / urea
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", hydrazine__CASES)
def test_hydrazine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
