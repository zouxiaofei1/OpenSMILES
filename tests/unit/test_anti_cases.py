# 合并自 4 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_anti_formic_branched_alkoxy.py: Reject formic-acid false parent for ring COOH; claim isopropoxy/isobutoxy on arenes.
test_anti_formyl_formate.py: Reject formaldehyde / methyl formate false parents for ring-only FG carbons.
test_anti_methanol.py: Reject methanol false parent for ring-carbon alcohols.
test_anti_methanone.py: Reject methanone false parent for ring-only ketone carbonyls.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_anti_formic_branched_alkoxy.py
# IUPAC: P-65.1.1.1 / P-63.2.2 / P-29.3
# Layer: L2,L3,L5
#
# Reject formic-acid false parent for ring COOH; claim isopropoxy/isobutoxy on arenes.
#
# Ring-only carboxyls (COOH carbon neighbors are ring/aryl only) must not collapse
# to open-chain formic acid. Branched outer alkoxy O–CHMe2 / O–CH2–CHMe2 must be
# claimable so simple benzoic gating passes.
# ==========================================================================
anti_formic_branched_alkoxy__POS_NO_FORMIC = [
    "IC1=C(SC(=C1S(=O)(=O)C(C)C)SC)C(=O)O",
    "CN1N=C(C=C1)C1=CC=C(C(=O)O)C=C1",
]

# Negative / regressions: true formic, open acids, methoxybenzoic
anti_formic_branched_alkoxy__NEG_CASES = [
    ("C(=O)O", "formic acid", "甲酸"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCCCCC(=O)O", "hexanoic acid", "己酸"),
    ("COc1ccc(C(=O)O)cc1", "4-methoxybenzoic acid", "4-甲氧基苯甲酸"),
]


@pytest.mark.parametrize("smiles", anti_formic_branched_alkoxy__POS_NO_FORMIC)
def test_ring_acid_not_formic(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "formic acid"
    assert "formic" not in en


@pytest.mark.parametrize("smiles,en,zh", anti_formic_branched_alkoxy__NEG_CASES)
def test_formic_and_open_acid_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_anti_formyl_formate.py
# IUPAC: P-66.6.1 / P-65.6 / P-44
# Layer: L2
#
# Reject formaldehyde / methyl formate false parents for ring-only FG carbons.
#
# Aldehyde or ester carbonyls that only attach to ring/aryl carbons must not
# collapse to open-chain C1 formaldehyde or methyl formate. True H2C=O and
# HCO2Me remain.
# ==========================================================================
anti_formyl_formate__POS_NO_FORMALDEHYDE = [
    "O=Cc1ccc(F)c(F)c1COc1ccccc1",
    "O=Cc1ccncc1",
]

anti_formyl_formate__POS_NO_METHYL_FORMATE = [
    "COC(=O)c1cccs1",  # methyl thiophene-2-carboxylate-like
]

anti_formyl_formate__NEG = [
    ("C=O", "formaldehyde", "甲醛"),
    ("COC=O", "methyl formate", "甲酸甲酯"),
]


@pytest.mark.parametrize("smiles", anti_formyl_formate__POS_NO_FORMALDEHYDE)
def test_ring_aldehyde_not_formaldehyde(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "formaldehyde"
    assert "formaldehyde" not in en


@pytest.mark.parametrize("smiles", anti_formyl_formate__POS_NO_METHYL_FORMATE)
def test_ring_ester_not_methyl_formate(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methyl formate"


@pytest.mark.parametrize("smiles,en,zh", anti_formyl_formate__NEG)
def test_formyl_formate_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_anti_methanol.py
# IUPAC: P-63.1.1 / P-44
# Layer: L2
#
# Reject methanol false parent for ring-carbon alcohols.
#
# OH on a ring carbon with no open-chain arm must not collapse to C1 methanol.
# Benzyl alcohol (exocyclic CH2OH) and true methanol must remain correct.
# ==========================================================================
anti_methanol__POS_NO_METHANOL = [
    "Nc1ccc2c(c1)CC(O)C2",  # aminoindanol-like
    "CC1CC(O)CC(C)(C)C1",  # substituted cyclohexanol not simple
]

anti_methanol__NEG_CASES = [
    ("CO", "methanol", "甲醇"),
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
]


@pytest.mark.parametrize("smiles", anti_methanol__POS_NO_METHANOL)
def test_ring_alcohol_not_methanol(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methanol"
    assert not en.endswith("methanol") or "phenyl" in en


@pytest.mark.parametrize("smiles,en,zh", anti_methanol__NEG_CASES)
def test_alcohol_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_anti_methanone.py
# IUPAC: P-64.2.1 / P-44
# Layer: L2
#
# Reject methanone false parent for ring-only ketone carbonyls.
#
# When the ketone carbon has no open-chain carbon arms (only ring/aryl C
# neighbors), open-chain kind=ketone must not collapse to C1 methanone.
# ==========================================================================
anti_methanone__POS_NO_METHANONE = [
    "O[C@]1(C(=CC(C1)=O)C1=CC=CC=C1)C1=CC=CC=C1",
    "C(CCC)C=1OC2=C(C(C1)=O)C=C(C=C2)Cl",
]

anti_methanone__NEG_CASES = [
    ("CCC(=O)CC", "pentan-3-one", "戊-3-酮"),
    ("CC(=O)C(C)C", "3-methylbutan-2-one", None),
]


@pytest.mark.parametrize("smiles", anti_methanone__POS_NO_METHANONE)
def test_ring_ketone_not_methanone(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methanone"
    assert "methanone" not in en


@pytest.mark.parametrize("smiles,en,zh", anti_methanone__NEG_CASES)
def test_open_ketone_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
