# 合并自 9 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_cycloalkane.py: Simple unsubstituted monocycloalkanes C3–C10 (P-22.1.1 / P-31).
test_mono_cycloalkene.py: Unsubstituted monocyclic monoalkenes: cyclohexene etc.
test_mono_cycloalcohol.py: Unsubstituted monocyclic monoalcohols (cycloalkanols).
test_mono_cycloamine.py: Unsubstituted monocyclic primary monoamines: cyclohexanamine etc.
test_mono_cycloketone.py: Unsubstituted monocyclic monoketones: cyclohexanone etc.
test_cycloalkyl.py: Monocycloalkyl substituents (unsubstituted C3–C8) on cycloalkane parents.
test_cycloalkyne_aryne.py: Endocyclic alkyne carbocycle (aryne): cyclohexa-1,3-dien-5-yne.
test_cyclopolyene.py: Unsubstituted monocyclic polyenes (cyclopolyene): cyclohexadiene etc.
test_cycloalkanediol_dione.py: Saturated monocyclic diols and diones (cycloalkanediol / cycloalkanedione).
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mono_cycloalkane.py
# IUPAC: P-22.1.1
# Layer: L1,L2,L5
#
# Simple unsubstituted monocycloalkanes C3–C10 (P-22.1.1 / P-31).
#
# cyclo + alkane stem; Chinese 环 + 烷. No side chains, no unsaturation.
# ==========================================================================
mono_cycloalkane__CASES = [
    ("C1CCC1", "cyclobutane", "环丁烷"),
    ("C1CCCC1", "cyclopentane", "环戊烷"),
    ("C1CCCCCC1", "cycloheptane", "环庚烷"),
    ("C1CCCCCCC1", "cyclooctane", "环辛烷"),
    ("C1CCCCCCCC1", "cyclononane", "环壬烷"),
    ("C1CCCCCCCCC1", "cyclodecane", "环癸烷"),
    ("CCCC", "butane", "丁烷"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cycloalkane__CASES)
def test_mono_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methylcyclohexane_not_cyclohexane() -> None:
    """Substituted ring is out of scope; must not be named cyclohexane."""
    r = SMILESNNamer().name("CC1CCCCC1")
    assert r.success
    assert normalize_en(r.en) != "cyclohexane"


# ==========================================================================
# 合并自 test_mono_cycloalkene.py
# IUPAC: P-31.1 / P-22.1.1
# Layer: L2,L4,L5
#
# Unsubstituted monocyclic monoalkenes: cyclohexene etc.
#
# Parent = monocarbocycle with one endocyclic C=C; no substituents.
# Unsubstituted: omit locant (cyclohexene not cyclohex-1-ene).
# ==========================================================================
mono_cycloalkene__CASES = [
    ("C1=CCC1", "cyclobutene", "环丁烯"),
    ("C1=CCCCCC1", "cycloheptene", "环庚烯"),
    ("NC1CCCCC1", "cyclohexanamine", "环己胺"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cycloalkene__CASES)
def test_mono_cycloalkene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_ene_name() -> None:
    r = SMILESNNamer().name("C1=CCCCC1")
    assert "dec" not in normalize_en(r.en)
    assert normalize_en(r.en) != "hexene"


# ==========================================================================
# 合并自 test_mono_cycloalcohol.py
# IUPAC: P-63.1.1 / P-22.1.1
# Layer: L2,L4,L5
#
# Unsubstituted monocyclic monoalcohols (cycloalkanols).
#
# Parent = saturated monocarbocycle with one ring-carbon OH.
# Unsubstituted: omit locant (cyclohexanol not cyclohexan-1-ol).
# ==========================================================================
mono_cycloalcohol__CASES = [
    # positive: unsubstituted monocycloalkanols C3–C7
    ("OC1CC1", "cyclopropanol", "环丙醇"),
    ("OC1CCC1", "cyclobutanol", "环丁醇"),
    ("OC1CCCC1", "cyclopentanol", "环戊醇"),
    ("OC1CCCCCC1", "cycloheptanol", "环庚醇"),
    ("CCCO", "propan-1-ol", "丙-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cycloalcohol__CASES)
def test_mono_cycloalcohol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cyclohexanol_not_chain_alcohol() -> None:
    """Ring OH must not be named as acyclic alcohol (hexanol / nonan-*-ol)."""
    r = SMILESNNamer().name("OC1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexanol"
    assert "hexanol" != en or en.startswith("cyclo")
    assert "nonan" not in en


# ==========================================================================
# 合并自 test_mono_cycloamine.py
# IUPAC: P-62.2.1 / P-22.1.1
# Layer: L2,L4,L5
#
# Unsubstituted monocyclic primary monoamines: cyclohexanamine etc.
#
# Parent = saturated monocarbocycle with one ring primary amine.
# Unsubstituted: omit locant (cyclohexanamine not cyclohexan-1-amine).
# ==========================================================================
mono_cycloamine__CASES = [
    ("NC1CCCC1", "cyclopentanamine", "环戊胺"),
    ("NC1CCC1", "cyclobutanamine", "环丁胺"),
    ("NC1CC1", "cyclopropanamine", "环丙胺"),
    ("NC1CCCCCC1", "cycloheptanamine", "环庚胺"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cycloamine__CASES)
def test_mono_cycloamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_amine_name() -> None:
    r = SMILESNNamer().name("NC1CCCCC1")
    assert "nonan" not in normalize_en(r.en)
    assert normalize_en(r.en) != "hexanamine"


# ==========================================================================
# 合并自 test_mono_cycloketone.py
# IUPAC: P-64.2.1 / P-22.1.1
# Layer: L2,L4,L5
#
# Unsubstituted monocyclic monoketones: cyclohexanone etc.
#
# Parent = saturated monocarbocycle with one ring carbonyl.
# Unsubstituted: omit locant (cyclohexanone not cyclohexan-1-one).
# ==========================================================================
mono_cycloketone__CASES = [
    ("CC(=O)CC", "butan-2-one", "丁-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_cycloketone__CASES)
def test_mono_cycloketone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_ketone_name() -> None:
    r = SMILESNNamer().name("O=C1CCCCC1")
    assert normalize_en(r.en) != "hexanone"
    assert "nonan" not in normalize_en(r.en)


# ==========================================================================
# 合并自 test_cycloalkyl.py
# IUPAC: P-29.6 / P-22.1.1 / P-14.3.4
# Layer: L2,L3,L5
#
# Monocycloalkyl substituents (unsubstituted C3–C8) on cycloalkane parents.
#
# P-29.6: cycloalkyl prefixes (cyclohexyl, cyclopentyl, …).
# Multi-ring molecules pick one sat carbocycle as parent; the other is a side.
# ==========================================================================
cycloalkyl__CASES = [
    # P0: bi(cyclohexane)
    ("C1CCC(CC1)C1CCCCC1", "cyclohexylcyclohexane", "环己基环己烷"),
    # P1: 1-cyclohexylethyl
    (
        "C1(CCCCC1)C(C)C1CCCCC1",
        "(1-cyclohexylethyl)cyclohexane",
        "(1-环己基乙基)环己烷",
    ),
    # generality: cyclopentyl
    ("C1CCC(CC1)C1CCCC1", "cyclopentylcyclohexane", "环戊基环己烷"),
    # regressions
    ("C1CCC(CC1)C", "methylcyclohexane", "甲基环己烷"),
    ("CC(C)(C)C1CCCCC1", "tert-butylcyclohexane", "叔丁基环己烷"),
    (
        "C1CC(C(C)(CC)C)CCC1",
        "(2-methylbutan-2-yl)cyclohexane",
        "(2-甲基丁-2-基)环己烷",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkyl__CASES)
def test_cycloalkyl(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_cycloalkyne_aryne.py
# IUPAC: P-31.1 / P-31.2
# Layer: L2,L4,L5
#
# Endocyclic alkyne carbocycle (aryne): cyclohexa-1,3-dien-5-yne.
#
# 环内三键使该环无法构成 mancude 芳香体系：RDKit 把苯炔 c1ccccc#1 的 sp 碳按 6π 一并标为
# 芳香，单环全芳香又匹配不到保留名时母体选择落空（no_assemblable_candidate）。环内三键
# 须按环烯炔表达，不得当保留芳名或未注册芳环丢弃。苯与纯烯环不得回归。
# ==========================================================================
cycloalkyne_aryne__CASES = [
    # 苯炔（C6H4）：三种等价 Kekulé 写法须同名
    ("C1=CC=CC#C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    ("C1=CC#CC=C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    ("C1#CC=CC=C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    # 带取代基的芳炔：环内三键走 carbocycle 路径，取代基照常表达
    ("CC1=CC=CC#C1", "1-methylcyclohexa-1,3-dien-5-yne", "1-甲基环己-1,3-二烯-5-炔"),
    # 非芳香环炔（RDKit 本就不标芳香）：不得受影响
    ("C1#CC=CC1", "cyclopent-1-en-3-yne", "环戊-1-烯-3-炔"),
    # 负例：真芳香环与纯烯环不得被当成烯炔
    ("C1=CC=CC=C1", "benzene", "苯"),
    ("C1=CCC=CC1", "cyclohexa-1,4-diene", "环己-1,4-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkyne_aryne__CASES)
def test_cycloalkyne_aryne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles} failed: {(r.meta or {}).get('reason')}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CC=CC#C1", "C1=CC#CC=C1", "C1#CC=CC=C1"])
def test_aryne_not_named_benzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert "benzene" not in normalize_en(r.en)
    assert "苯" != normalize_zh(r.zh)


# ==========================================================================
# 合并自 test_cyclopolyene.py
# IUPAC: P-31.1 / P-22.1.3
# Layer: L2,L4,L5
#
# Unsubstituted monocyclic polyenes (cyclopolyene): cyclohexadiene etc.
#
# Parent = mono all-C carbocycle with ≥2 endocyclic non-aromatic C=C; no main FG.
# Ene locant set lowest (P-31.1); EN cyclohexa-1,3-diene; ZH 环己-1,3-二烯.
# Mono cycloalkene and open polyene must not regress.
# ==========================================================================
cyclopolyene__CASES = [
    # positive: ring dienes
    ("C1=CC=CCC1", "cyclohexa-1,3-diene", "环己-1,3-二烯"),
    ("C1=CCC=CC1", "cyclohexa-1,4-diene", "环己-1,4-二烯"),
    ("C1=CC=CC1", "cyclopenta-1,3-diene", "环戊-1,3-二烯"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", cyclopolyene__CASES)
def test_cyclopolyene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CC=CCC1", "C1=CCC=CC1", "C1=CC=CC1"])
def test_not_methane_or_alkane(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert en != "methane"
    assert "diene" in en
    assert "二烯" in (r.zh or "")


# ==========================================================================
# 合并自 test_cycloalkanediol_dione.py
# IUPAC: P-63.1.2 / P-64.2.1 / P-14.3.4
# Layer: L2,L4,L5
#
# Saturated monocyclic diols and diones (cycloalkanediol / cycloalkanedione).
#
# Parent = unfused C3–C10 sat carbocycle with exactly two ring-carbon OH
# or two ring ketone carbonyls. Lowest locant pair for the two FGs (always
# shown). Optional simple ring halo / claimable n-alkyl. Must not capture
# benzenediol, mono cycloalcohol/one, or open-chain diol/dione.
# ==========================================================================
cycloalkanediol_dione__CASES = [
    # positive: unsubstituted cycloalkanediols
    ("OC1CCCCC1O", "cyclohexane-1,2-diol", "环己烷-1,2-二醇"),
    ("OC1CC(O)CCC1", "cyclohexane-1,3-diol", "环己烷-1,3-二醇"),
    ("OC1CCC(O)CC1", "cyclohexane-1,4-diol", "环己烷-1,4-二醇"),
    ("OC1CCCC1O", "cyclopentane-1,2-diol", "环戊烷-1,2-二醇"),
    ("OC1CCC1O", "cyclobutane-1,2-diol", "环丁烷-1,2-二醇"),
    # positive: simple ring sub keeps FG pair + sub locant
    ("CC1CC(O)CC(O)C1", "5-methylcyclohexane-1,3-diol", "5-甲基环己烷-1,3-二醇"),
    ("ClC1CC(O)CC(O)C1", "5-chlorocyclohexane-1,3-diol", "5-氯环己烷-1,3-二醇"),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkanediol_dione__CASES)
def test_cycloalkanediol_dione(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cyclohexanediol_not_methane() -> None:
    """Ring 1,2-diol must not collapse to methane / open-chain hexanediol."""
    r = SMILESNNamer().name("OC1CCCCC1O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexane-1,2-diol"
    assert en != "methane"
    assert "benzene" not in en
