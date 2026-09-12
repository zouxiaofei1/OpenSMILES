# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_benzoic_acid.py: Retained parent name benzoic acid (benzene + one carboxyl on ring carbon).
test_aminobenzoic.py: Arene retained acid/aldehyde: ring amino/alkoxy/nitro as prefixes.
test_cycloalkanecarboxylic.py: Monocyclic cycloalkanecarboxylic acids (unsub or mono ring sub).
test_ring_polycarboxylic.py: 环烷/芳环上多个 exocyclic COOH 的系统命名。
test_hetero5_carboxylic.py: Hetero5 carboxylic acid scope (P-65.1.1) — negative guard only.
test_pyridinecarboxylic_zh.py: Pyridinecarboxylic Chinese suffix scope (P-65.1.1) — negative guard only.
test_sat_hetero_carboxylic.py: Sat-monohetero carboxylic acid scope (P-65.1.1) — negative guard only.
test_indole_naph_carboxylic.py: Retained fused-ring carboxylic parents: indolecarboxylic / naphthalenecarboxylic.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_benzoic_acid.py
# IUPAC: P-65.1.1.1
# Layer: L2,L3,L4,L5
#
# Retained parent name benzoic acid (benzene + one carboxyl on ring carbon).
#
# Parent = benzene with principal characteristic COOH (benzoic acid).
# Allow 0–2 extra ring halo / methyl / phenolic hydroxy; COOH attach = locant 1.
# ==========================================================================
benzoic_acid__CASES = [
    # positive: unsubstituted retained parent
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    # positive: monohalo / monomethyl
    ("OC(=O)c1ccc(Cl)cc1", "4-chlorobenzoic acid", "4-氯苯甲酸"),
    ("Cc1cccc(C(=O)O)c1", "3-methylbenzoic acid", None),
    # positive: phenolic hydroxy as prefix
    ("OC(=O)c1ccc(O)cc1", "4-hydroxybenzoic acid", "4-羟基苯甲酸"),
    ("O=C(O)c1cccc(O)c1", "3-hydroxybenzoic acid", None),
]


@pytest.mark.parametrize("smiles,en,zh", benzoic_acid__CASES)
def test_benzoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzoic_not_chain_acid() -> None:
    """Aromatic carboxyl must not expand ring into chain acid parent."""
    r = SMILESNNamer().name("OC(=O)c1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzoic acid"
    assert "heptanoic" not in en
    assert "hexanoic" not in en


# ==========================================================================
# 合并自 test_aminobenzoic.py
# IUPAC: P-65.1.1.1
# Layer: L2,L3,L5
#
# Arene retained acid/aldehyde: ring amino/alkoxy/nitro as prefixes.
#
# Benzoic acid (and benzaldehyde) keep COOH/CHO as principal FG; simple
# ring amino, methoxy/ethoxy, and nitro are allowed prefixes alongside
# existing halo/methyl/hydroxy (sub cap ≤3). Aligns with phenol/benzene.
# ==========================================================================
aminobenzoic__CASES = [
    # positive: aminobenzoic
    ("Nc1ccccc1C(=O)O", "2-aminobenzoic acid", None),
    ("Nc1ccc(C(=O)O)cc1", "4-aminobenzoic acid", "4-氨基苯甲酸"),
    # positive: amino + methoxy + dimethyl
    ("NC1=CC(=C(C(=O)O)C=C1)OC", "4-amino-2-methoxybenzoic acid", "4-氨基-2-甲氧基苯甲酸"),
    ("NC1=C(C=C(C(=O)O)C=C1C)C", "4-amino-3,5-dimethylbenzoic acid", None),
    ("O=[N+]([O-])c1ccc(C(=O)O)cc1", "4-nitrobenzoic acid", "4-硝基苯甲酸"),
    # optional benzaldehyde
    ("O=Cc1ccc(OC)cc1", "4-methoxybenzaldehyde", "4-甲氧基苯甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", aminobenzoic__CASES)
def test_aminobenzoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_aminobenzoic_not_chain_acid() -> None:
    """Ring carboxyl + amino must not expand to aminoalkanoic acid."""
    r = SMILESNNamer().name("Nc1ccccc1C(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-aminobenzoic acid"
    assert "heptanoic" not in en
    assert "hexanoic" not in en


# ==========================================================================
# 合并自 test_cycloalkanecarboxylic.py
# IUPAC: P-65.1.1
# Layer: L2,L4,L5
#
# Monocyclic cycloalkanecarboxylic acids (unsub or mono ring sub).
#
# Parent = saturated monocarbocycle with one exocyclic COOH on a ring carbon.
# Allow 0–1 ring simple sub: monohalo (F/Cl/Br/I) or mono C1–C4 n-alkyl.
# COOH attachment = locant 1; ring subs get lowest compatible locants.
# IUPAC: cycloalkanecarboxylic acid / 环…烷甲酸.
# ==========================================================================
cycloalkanecarboxylic__CASES = [
    ("OC(=O)CCCCCC", "heptanoic acid", "庚酸"),
    ("CCCCCCC(=O)O", "heptanoic acid", "庚酸"),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkanecarboxylic__CASES)
def test_cycloalkanecarboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_ring_polycarboxylic.py
# IUPAC: P-65.2.2 (环/稠环上的多 exocyclic -COOH → -Xcarboxylic acid)
# Layer: L4,L5
#
# 环烷/芳环上多个 exocyclic COOH 的系统命名。
#
# Parent = carbocycle / benzene with multiplicity≥2 exocyclic COOH.
# 多羧酸用 -Xcarboxylic acid（非链式 Xanedioic 模板），必须带全部位次:
#   cyclohexane-1,2-dicarboxylic acid / 环己烷-1,2-二羧酸
#   benzene-1,4-dicarboxylic acid        / 苯-1,4-二羧酸
# 单羧酸（cyclohexanecarboxylic acid / benzoic acid）不得回归。
# ==========================================================================
ring_polycarboxylic__CASES = [
    # 环己烷-1,2-二甲酸（用户示例，ortho）
    ("C1CCC(C(=O)O)C(C(=O)O)C1",
     "cyclohexane-1,2-dicarboxylic acid", "环己烷-1,2-二羧酸"),
    # 环己烷-1,4-二甲酸（para）
    ("O=C(O)C1CCC(C(=O)O)CC1",
     "cyclohexane-1,4-dicarboxylic acid", "环己烷-1,4-二羧酸"),
    # 苯-1,2-二甲酸（phthalic）
    ("C1=CC=C(C(=O)O)C(=C1)C(=O)O",
     "benzene-1,2-dicarboxylic acid", "苯-1,2-二羧酸"),
    # 苯-1,4-二甲酸（对苯二甲酸系统名）
    ("O=C(O)c1ccc(C(=O)O)cc1",
     "benzene-1,4-dicarboxylic acid", "苯-1,4-二羧酸"),
    # 苯-1,3,5-三甲酸（trimesic 系统名）
    ("O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1",
     "benzene-1,3,5-tricarboxylic acid", "苯-1,3,5-三羧酸"),
]

# 单羧酸对照：不得回归
ring_polycarboxylic__NEG = [
    ("C1CCC(C(=O)O)CC1", "cyclohexanecarboxylic acid", "环己烷羧酸"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", ring_polycarboxylic__CASES)
def test_ring_polycarboxylic(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


@pytest.mark.parametrize("smiles,en,zh", ring_polycarboxylic__NEG)
def test_monocarboxylic_no_regress(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


# ==========================================================================
# 合并自 test_hetero5_carboxylic.py
# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
#
# Hetero5 carboxylic acid scope (P-65.1.1) — negative guard only.
#
# Positive hetero5-acid cases are not covered in this file. The retained case
# asserts benzoic acid is not renamed as a hetero5 carboxylic acid.
# ==========================================================================
hetero5_carboxylic__CASES = [
    # negative: must not steal existing correct parents
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", hetero5_carboxylic__CASES)
def test_hetero5_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyridinecarboxylic_zh.py
# IUPAC: P-65.1.1
# Layer: L5
#
# Pyridinecarboxylic Chinese suffix scope (P-65.1.1) — negative guard only.
#
# Positive 吡啶-n-甲酸 cases are not covered in this file. The retained case
# asserts benzoic acid keeps 苯甲酸 (no regression to a pyridine suffix).
# ==========================================================================
pyridinecarboxylic_zh__CASES = [
    # negative: carbocyclic carboxylic acids keep correct 甲酸 (no regression)
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridinecarboxylic_zh__CASES)
def test_pyridinecarboxylic_zh(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_sat_hetero_carboxylic.py
# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
#
# Sat-monohetero carboxylic acid scope (P-65.1.1) — negative guard only.
#
# Positive sat-hetero acid cases are not covered in this file. The retained cases
# assert benzoic acid and hexanoic acid are not renamed as sat-hetero acids.
# ==========================================================================
sat_hetero_carboxylic__CASES = [
    ("O=C(O)CCCCC", "hexanoic acid", "己酸"),
]


@pytest.mark.parametrize("smiles,en,zh", sat_hetero_carboxylic__CASES)
def test_sat_hetero_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_indole_naph_carboxylic.py
# IUPAC: P-65.1.1 / P-22.2.1 / P-25
# Layer: L2,L4,L5
#
# Retained fused-ring carboxylic parents: indolecarboxylic / naphthalenecarboxylic.
#
# IUPAC P-65.1.1 / P-22.2.1 / P-25: mono COOH on retained indole or naphthalene
# core yields 1H-indole-n-carboxylic acid / naphthalene-n-carboxylic acid
# (Chinese: 1H-吲哚-n-甲酸 / 萘-n-甲酸). Must not collapse to open-chain formic acid.
# ==========================================================================
indole_naph_carboxylic__CASES = [
    ("O=C(O)[C@@H](O)CO", "(2S)-2,3-dihydroxypropanoic acid", "(2S)-2,3-二羟基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", indole_naph_carboxylic__CASES)
def test_indole_naph_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
