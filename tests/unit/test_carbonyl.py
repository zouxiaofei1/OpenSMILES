# 合并自 13 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_aldehyde.py: Simple acyclic monoaldehydes (alkanals): retained C1/C2 + systematic C3+.
test_mono_ketone.py: Simple acyclic monoketones (alkanones): stem + locant + one/酮.
test_alkenal.py: Open-chain monounsaturated monoaldehydes (alkenals).
test_alkenone.py: Open-chain monounsaturated monoketones (alkenones / alkynones).
test_alkanedione.py: Unsubstituted open-chain saturated alkanediones (exactly two ketones).
test_polyalkenal.py:
test_benzaldehyde.py: Retained parent name benzaldehyde (benzene + one formyl on ring carbon).
test_acetophenone.py: Acetophenone scope (P-64.1.1) — negative guard only.
test_anthraquinone.py: 9,10-Anthraquinone retained parent (linear anthracene + meso dione).
test_benzoquinone.py: 1,4-Benzoquinone PIN as cyclohexa-2,5-diene-1,4-dione (P-64.2).
test_ortho_benzoquinone.py: 1,2-Benzoquinone PIN as cyclohexa-3,5-diene-1,2-dione (P-64.2 Round C).
test_chromenone.py: Retained parent chromen-2-one (coumarin lactone; 2H-1-benzopyran-2-one style).
test_ring_carbaldehyde_formyl.py: 环上外环 -CHO：作主官能团时 -carbaldehyde（单/多），降级为取代基前缀时 formyl。
"""
from __future__ import annotations

import pytest

from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mono_aldehyde.py
# IUPAC: P-66.6.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic monoaldehydes (alkanals): retained C1/C2 + systematic C3+.
#
# Aldehyde carbonyl carbon has exactly one carbon neighbor and is not carboxyl.
# Suffix -al / 醛 has no locant (aldehyde carbon is always position 1).
# ==========================================================================
mono_aldehyde__CASES = [
    ("CCC=O", "propanal", "丙醛"),
    ("CCCCCC=O", "hexanal", "己醛"),
    ("CC(C)C=O", "2-methylpropanal", "2-甲基丙醛"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_aldehyde__CASES)
def test_mono_aldehyde(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_mono_ketone.py
# IUPAC: P-64.2.1
# Layer: L1,L2,L4,L5
#
# Simple acyclic monoketones (alkanones): stem + locant + one/酮.
#
# Carbonyl carbon has exactly two carbon neighbors and is not carboxyl.
# ==========================================================================
mono_ketone__CASES = [
    ("CCC(CC)=O", "pentan-3-one", "戊-3-酮"),
    ("CCCCC(=O)CC", "heptan-3-one", "庚-3-酮"),
    ("CCCCCC(=O)C", "heptan-2-one", "庚-2-酮"),
    ("CCCCCCCCC(C)=O", "decan-2-one", None),
]


@pytest.mark.parametrize("smiles,en,zh", mono_ketone__CASES)
def test_mono_ketone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenal.py
# IUPAC: P-66.6.1 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated monoaldehydes (alkenals).
#
# Aldehyde is the principal characteristic group (CHO = locant 1); one non-aromatic
# C=C is expressed as -n-enal / -n-烯醛 with the lower double-bond carbon locant.
# No (E)/(Z) stereodescriptors this cycle.
# ==========================================================================
alkenal__CASES = [
    # positive: acyclic mono-alkenals (no E/Z)
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("C=CC=O", "prop-2-enal", "丙-2-烯醛"),
    ("CCC=CC=O", "pent-2-enal", "戊-2-烯醛"),
    ("CCCC=CC=O", "hex-2-enal", "己-2-烯醛"),
    ("C=CCC=O", "but-3-enal", "丁-3-烯醛"),
    ("CCCC=O", "butanal", "丁醛"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenal__CASES)
def test_alkenal(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenone.py
# IUPAC: P-64.2.1 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated monoketones (alkenones / alkynones).
#
# Ketone is the principal characteristic group; one non-aromatic C=C or C≡C is
# expressed as -n-en-m-one / -n-烯-m-酮 or -n-yn-m-one / -n-炔-m-酮. Ketone locant
# is minimized first (P-64.2.1); when tied, unsaturation takes the lower set
# (P-31.1). E/Z via BondStereo when defined.
# ==========================================================================
alkenone__CASES = [
    # positive: acyclic mono-alkenones
    ("C=CC(C)=O", "but-3-en-2-one", "丁-3-烯-2-酮"),
    ("CC=CC(C)=O", "pent-3-en-2-one", "戊-3-烯-2-酮"),
    ("C=CCC(C)=O", "pent-4-en-2-one", "戊-4-烯-2-酮"),
    ("CCC=CC(C)=O", "hex-3-en-2-one", "己-3-烯-2-酮"),
    # positive: alkynone
    ("C#CC(C)=O", "but-3-yn-2-one", "丁-3-炔-2-酮"),
    # E/Z mono-alkenone
    (r"C/C=C/C(C)=O", "(3E)-pent-3-en-2-one", "(3E)-戊-3-烯-2-酮"),
    # negatives: saturated ketone / ring / aromatic must not become alkenone
    ("CCC(C)=O", "butan-2-one", "丁-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenone__CASES)
def test_alkenone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanedione.py
# IUPAC: P-64.2.1
# Layer: L2,L4,L5
#
# Unsubstituted open-chain saturated alkanediones (exactly two ketones).
#
# Parent chain through both ketone carbons; name alkane-{a},{b}-dione /
# {首字}-{a},{b}-二酮 (mono-ketone Chinese style).
# ==========================================================================
alkanedione__CASES = [
    ("CC(=O)C(C)=O", "butane-2,3-dione", "丁-2,3-二酮"),
    ("CC(=O)CCC(=O)C", "hexane-2,5-dione", "己-2,5-二酮"),
    ("CC(=O)CCCCC(=O)C", "octane-2,7-dione", "辛-2,7-二酮"),
    ("CCC(=O)CC(=O)CC", "heptane-3,5-dione", "庚-3,5-二酮"),
    # positive: trione (multiplicity-generic count suffix)
    ("CC(=O)CC(=O)CC(=O)C", "heptane-2,4,6-trione", "庚-2,4,6-三酮"),
    ("CCCC(=O)C", "pentan-2-one", "戊-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanedione__CASES)
def test_alkanedione(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_polyalkenal.py
# IUPAC: P-66.6.1 / P-31.1 / P-93.4
# Layer: L5
# ==========================================================================
polyalkenal__CASES = [
    ("C=CC=CC=O", "penta-2,4-dienal", "戊-2,4-二烯醛"),
    ("CC=CC=CC=O", "hexa-2,4-dienal", "己-2,4-二烯醛"),
    ("CC(C)=CCC/C(C)=C/C=O", "(2E)-3,7-dimethylocta-2,6-dienal", "(2E)-3,7-二甲基辛-2,6-二烯醛"),
    # negatives
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("CCCC=O", "butanal", "丁醛"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
]

@pytest.mark.parametrize("smiles,en,zh", polyalkenal__CASES)
def test_polyalkenal(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzaldehyde.py
# IUPAC: P-66.6.1
# Layer: L2,L4,L5
#
# Retained parent name benzaldehyde (benzene + one formyl on ring carbon).
#
# Parent = benzene with principal characteristic CHO (benzaldehyde).
# Allow 0–2 extra ring halo / methyl / phenolic hydroxy; CHO attach = locant 1.
# ==========================================================================
benzaldehyde__CASES = [

    # positive: phenolic hydroxy as prefix
    ("O=Cc1ccc(O)cc1", "4-hydroxybenzaldehyde", None),
]


@pytest.mark.parametrize("smiles,en,zh", benzaldehyde__CASES)
def test_benzaldehyde(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzaldehyde_not_chain_al() -> None:
    """Aromatic formyl must not expand ring into chain aldehyde parent."""
    r = SMILESNNamer().name("O=Cc1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzaldehyde"
    assert "heptanal" not in en
    assert "hexanal" not in en


# ==========================================================================
# 合并自 test_acetophenone.py
# IUPAC: P-64.1.1
# Layer: L2,L4,L5
#
# Acetophenone scope (P-64.1.1) — negative guard only.
#
# Positive acetophenone cases are not covered in this file. The retained cases
# assert benzaldehyde / propan-2-one / benzene / benzoic acid / ethanol are not
# named as acetophenone.
# ==========================================================================
acetophenone__CASES = [
    # negative: benzaldehyde, acetone, benzene, benzoic acid, ethanol
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", acetophenone__CASES)
def test_acetophenone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_anthraquinone.py
# IUPAC: P-25 / P-64
# Layer: L2,L4,L5
#
# 9,10-Anthraquinone retained parent (linear anthracene + meso dione).
#
# PIN-style general name: 9,10-anthraquinone / 蒽醌 with ring methyl prefixes.
# ==========================================================================
anthraquinone__CASES = [
    # open-chain alkanedione not captured as anthraquinone
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
]


@pytest.mark.parametrize("smiles,en,zh", anthraquinone__CASES)
def test_anthraquinone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzoquinone.py
# IUPAC: P-64.2
# Layer: L2,L3
#
# 1,4-Benzoquinone PIN as cyclohexa-2,5-diene-1,4-dione (P-64.2).
#
# Single C6 carbocycle + exactly two para ring ketones + two endocyclic
# double bonds. PIN is systematic (not retained 1,4-benzoquinone).
#
# Round B: ring n-alkyl C1–C12 + halo/OH/alkoxy C1–C2 mixed substitution.
# ==========================================================================



# ==========================================================================
# 合并自 test_ortho_benzoquinone.py
# IUPAC: P-64.2
# Layer: L2,L4,L5
#
# 1,2-Benzoquinone PIN as cyclohexa-3,5-diene-1,2-dione (P-64.2 Round C).
#
# Single C6 carbocycle + exactly two ortho ring ketones + two endocyclic C=C.
# PIN is systematic (not retained o-benzoquinone). First-cut: halo / OH / Me
# or n-alkyl C1–C12 / alkoxy C1–C2. No naphthoquinone this round.
# ==========================================================================
ortho_benzoquinone__CASES = [
    # unsubstituted parent (was 1,4-BQ negative methanone)
    (
        "O=C1C=CC=CC1=O",
        "cyclohexa-3,5-diene-1,2-dione",
        "环己-3,5-二烯-1,2-二酮",
    ),
    # gold: 4-hydroxy-5-methyl
    (
        "CC1=CC(=O)C(=O)C=C1O",
        "4-hydroxy-5-methylcyclohexa-3,5-diene-1,2-dione",
        "4-羟基-5-甲基环己-3,5-二烯-1,2-二酮",
    ),
    # simple methyl
    (
        "CC1=CC=CC(=O)C1=O",
        "3-methylcyclohexa-3,5-diene-1,2-dione",
        "3-甲基环己-3,5-二烯-1,2-二酮",
    ),
    # 4,5-dichloro
    (
        "ClC1=CC(=O)C(=O)C=C1Cl",
        "4,5-dichlorocyclohexa-3,5-diene-1,2-dione",
        "4,5-二氯环己-3,5-二烯-1,2-二酮",
    ),
    # negatives: 1,4-BQ must stay; anthraquinone / acetophenone unchanged
    (
        "O=C1C=CC(=O)C=C1",
        "cyclohexa-2,5-diene-1,4-dione",
        "环己-2,5-二烯-1,4-二酮",
    ),
    ("CC(=O)c1ccccc1", "acetophenone", "苯乙酮"),
    ("O=C1c2ccccc2C(=O)c2ccccc12", "9,10-anthraquinone", "蒽醌"),
]


# ==========================================================================
# 合并自 test_chromenone.py
# IUPAC: P-25 / P-22.2.1 / P-65.6.3
# Layer: L2,L3,L4,L5
#
# Retained parent chromen-2-one (coumarin lactone; 2H-1-benzopyran-2-one style).
#
# Fused 6+6: benzene + α-pyrone lactone (9C + ring O + exocyclic =O at C2).
# O=1, carbonyl C=2; EN stem chromen-2-one; ZH 香豆素. Simple ring prefixes
# (halo / Me / alkoxy / OH / amino / nitro / prenyl / phenyl) only.
# ==========================================================================
chromenone__CASES = [
    # positive: unsubstituted retained parent
    ("O=c1ccc2ccccc2o1", "chromen-2-one", "香豆素"),
    # positive: simple ring prefixes
    ("COc1ccc2ccc(=O)oc2c1", "7-methoxychromen-2-one", "7-甲氧基香豆素"),
    ("Cc1ccc2ccc(=O)oc2c1", "7-methylchromen-2-one", "7-甲基香豆素"),
    ("Clc1ccc2ccc(=O)oc2c1", "7-chlorochromen-2-one", "7-氯香豆素"),
    # positive: anchor (prenyl + methoxy)
    (
        "COc1cc2oc(=O)ccc2cc1CC=C(C)C",
        "7-methoxy-6-(3-methylbut-2-enyl)chromen-2-one",
        "7-甲氧基-6-(3-甲基丁-2-烯基)香豆素",
    ),
    # positive: optional phenyl (omit 2H- to match project stem)
    ("c1ccc(cc1)c1cc(=O)oc2ccccc12", "4-phenylchromen-2-one", "4-苯基香豆素"),
    # negative: must not become chromen-2-one
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
]


# ==========================================================================
# 合并自 test_ring_carbaldehyde_formyl.py
# IUPAC: P-66.6.1.1.3 (-CHO 连环/环系 → -carbaldehyde)、P-66.6.1.1(3) (-CHO 作取代基前缀 → formyl)
# Layer: L2,L4,L5 (+ L3 formyl 前缀 / anchored_table)
#
# 环上外环 -CHO：作主官能团时 -carbaldehyde（单/多），降级为取代基前缀时 formyl。
#
# 覆盖环骨架（吡啶/苯/环烷）、多基团 -dicarbaldehyde，及保留名
# benzaldehyde 不回归；formyl 前缀取代 oxomethyl（P-66.6.1.1(3)）。
# ==========================================================================
ring_carbaldehyde_formyl__CARBALDEHYDE = [
    # 吡啶环 + 甲氨基取代，醛 locant 由 N 起算定向（tiers-47927 黄金名）
    ("CNC1=C(C=O)C=CC=N1", "2-(methylamino)pyridine-3-carbaldehyde", "2-(甲氨基)吡啶-3-甲醛"),
    ("O=CC1=CC=CC=N1", "pyridine-2-carbaldehyde", "吡啶-2-甲醛"),
    ("O=CC1=NC(C=O)=CC=C1", "pyridine-2,6-dicarbaldehyde", "吡啶-2,6-二甲醛"),
    # 苯环多醛 → 系统名 dicarbaldehyde（非保留 phthalaldehyde）
    ("C1C(C=O)=C(C=O)C=CC=1", "benzene-1,2-dicarbaldehyde", "苯-1,2-二甲醛"),
    ("O=Cc1ccc(C=O)cc1", "benzene-1,4-dicarbaldehyde", "苯-1,4-二甲醛"),
    # 单环环烷烃：醛基位次 1 隐含省略
    ("O=CC1CCCCC1", "cyclohexanecarbaldehyde", "环己烷甲醛"),
    # 哒嗪-3-甲醛（两个相邻等价环 N，方向不能随 SMILES 书写顺序变，
    # P-14.4(c)：后缀醛占低位次 3，CF3 前缀只得 6 —— 同分子两种写法必须同名）
    ("FC(C1=CC=C(N=N1)C=O)(F)F", "6-(trifluoromethyl)pyridazine-3-carbaldehyde", "6-(三氟甲基)哒嗪-3-甲醛"),
    ("C1(C(F)(F)F)N=NC(C=O)=CC=1", "6-(trifluoromethyl)pyridazine-3-carbaldehyde", "6-(三氟甲基)哒嗪-3-甲醛"),
]

ring_carbaldehyde_formyl__FORMYL_PREFIX = [
    ("OC(=O)c1ccc(cc1)C=O", "4-formylbenzoic acid", "4-甲酰苯甲酸"),
    ("OC(=O)c1ccccc1C=O", "2-formylbenzoic acid", "2-甲酰苯甲酸"),
    ("NC(=O)c1ccc(C=O)cc1", "4-formylbenzamide", "4-甲酰苯甲酰胺"),
]

# 保留名 / 开链醛不回归（P-66.6.1.2 苯甲醛系保留名；不误伤链醛）。
ring_carbaldehyde_formyl__REGRESSION = [
    ("Cc1ccc(C=O)cc1", "4-methylbenzaldehyde", "4-甲基苯甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", ring_carbaldehyde_formyl__CARBALDEHYDE)
def test_ring_carbaldehyde(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", ring_carbaldehyde_formyl__FORMYL_PREFIX)
def test_formyl_prefix(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", ring_carbaldehyde_formyl__REGRESSION)
def test_regression(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
