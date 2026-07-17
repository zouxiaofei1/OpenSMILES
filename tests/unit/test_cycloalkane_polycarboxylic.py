# IUPAC: P-65.1.1/P-65.1.2
# Layer: L2,L4,L5
"""Neutral, unfused saturated cycloalkane di- and tricarboxylic acids."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer2.cyclo_polycarboxylic import _try_cycloalkane_polycarboxylic
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.numbering import number
from namepredict.layer4.locants.plan import NumberingPlan
from namepredict.namer import SMILESNNamer


CASES = [
    ("O=C(O)C1CCCCC1C(=O)O", "cyclohexane-1,2-dicarboxylic acid", "环己烷-1,2-二甲酸"),
    ("O=C(O)C1CCCC(C(=O)O)C1", "cyclohexane-1,3-dicarboxylic acid", "环己烷-1,3-二甲酸"),
    ("O=C(O)C1CCC(C(=O)O)CC1", "cyclohexane-1,4-dicarboxylic acid", "环己烷-1,4-二甲酸"),
    ("O=C(O)C1C(C(=O)O)C(C(=O)O)CCC1", "cyclohexane-1,2,3-tricarboxylic acid", "环己烷-1,2,3-三甲酸"),
    ("O=C(O)C1CC1C(=O)O", "cyclopropane-1,2-dicarboxylic acid", "环丙烷-1,2-二甲酸"),
    ("O=C(O)C1CCC1C(=O)O", "cyclobutane-1,2-dicarboxylic acid", "环丁烷-1,2-二甲酸"),
    ("O=C(O)C1CCCC1C(=O)O", "cyclopentane-1,2-dicarboxylic acid", "环戊烷-1,2-二甲酸"),
    ("O=C(O)C1CCCCCC1C(=O)O", "cycloheptane-1,2-dicarboxylic acid", "环庚烷-1,2-二甲酸"),
    ("O=C(O)C1CCCCCCC1C(=O)O", "cyclooctane-1,2-dicarboxylic acid", "环辛烷-1,2-二甲酸"),
    ("O=C(O)C1CCCCCCCC1C(=O)O", "cyclononane-1,2-dicarboxylic acid", "环壬烷-1,2-二甲酸"),
    ("O=C(O)C1CCCCCCCCC1C(=O)O", "cyclodecane-1,2-dicarboxylic acid", "环癸烷-1,2-二甲酸"),
    ("O=C(O)C1CC(C)CCC1C(=O)O", "4-methylcyclohexane-1,2-dicarboxylic acid", "4-甲基环己烷-1,2-二甲酸"),
    ("O=C(O)C1CC(Cl)CCC1C(=O)O", "4-chlorocyclohexane-1,2-dicarboxylic acid", "4-氯环己烷-1,2-二甲酸"),
]


def _numbered(smiles: str) -> dict:
    from namepredict.layer0.preprocessor import preprocess
    from namepredict.layer1.analyzer import analyze
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    parent = _try_cycloalkane_polycarboxylic(info)
    assert parent is not None
    return number(parent, extract_substituents(info, parent))


def test_cyclo_polyacid_uses_scaffold_numbering_plan() -> None:
    numbered = _numbered("O=C(O)C1CC(C)CCC1C(=O)O")
    plan = numbered["parent"].get("numbering")
    assert isinstance(plan, NumberingPlan)
    assert plan.scaffold_id == "cycloalkane_polycarboxylic"
    assert numbered["parent"]["stem_en"] == "cyclohexane"
    assert numbered["cooh_locants"] == (1, 2)
    assert numbered["substituents"][0]["locant"] == 4


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkane_polycarboxylic_acids(smiles: str, en: str, zh: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_cycloalkane_polyacid_is_invariant_to_smiles_order() -> None:
    expected = "cyclohexane-1,2-dicarboxylic acid"
    for smiles in ("O=C(O)C1CCCCC1C(=O)O", "O=C(O)C1CCCCC1C(O)=O"):
        result = SMILESNNamer().name(smiles)
        assert result.success
        assert normalize_en(result.en) == normalize_en(expected)


@pytest.mark.parametrize("smiles", [
    "O=C([O-])C1CCCCC1C(=O)O",
    "O=C(O)C1C(C(=O)O)C(C(=O)O)C(C(=O)O)CC1",
    "O=C(O)C1CC(C)C(C)CC1C(=O)O",
    "O=C(O)[C@H]1CCCCC1C(=O)O",
])
def test_cycloalkane_polyacid_exclusions_are_rejected(smiles: str) -> None:
    assert not SMILESNNamer().name(smiles).success


def test_benzene_polyacid_bypasses_cycloalkane_scope() -> None:
    result = SMILESNNamer().name("O=C(O)c1ccccc1C(=O)O")
    assert result.success
    assert normalize_en(result.en) == normalize_en("benzene-1,2-dicarboxylic acid")


@pytest.mark.parametrize(("smiles", "en", "zh"), [
    ("O=C(O)[C@H]1CCC[C@H](C(=O)O)C1", "trans-cyclohexane-1,3-dicarboxylic acid", "反-环己烷-1,3-二甲酸"),
    ("O=C(O)[C@H]1CCC[C@@H](C(=O)O)C1", "cis-cyclohexane-1,3-dicarboxylic acid", "顺-环己烷-1,3-二甲酸"),
])
def test_cyclo_polyacid_two_site_relative_stereo(smiles: str, en: str, zh: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_cyclo_polyacid_three_site_relative_stereo_uses_plan_locants() -> None:
    result = SMILESNNamer().name("O=C(O)[C@H]1[C@@H](C(=O)O)C[C@H](C(=O)O)CC1")
    assert result.success
    assert normalize_en(result.en) == normalize_en("1r,2c,4t-cyclohexane-1,2,4-tricarboxylic acid")
    assert normalize_zh(result.zh) == normalize_zh("1r,2c,4t-环己烷-1,2,4-三甲酸")


def test_cyclo_polyacid_relative_name_is_enantiomer_invariant() -> None:
    names = [SMILESNNamer().name(s).en for s in (
        "O=C(O)[C@H]1CCC[C@@H](C(=O)O)C1",
        "O=C(O)[C@@H]1CCC[C@H](C(=O)O)C1",
    )]
    assert names == ["cis-cyclohexane-1,3-dicarboxylic acid"] * 2


@pytest.mark.parametrize("smiles", [
    "O=C(O)[C@H]1CCCCC1C(=O)O",
    "O=C(O)[C@H]1CCC[C@H](C(=O)O)C1Cl",
])
def test_cyclo_polyacid_incomplete_or_out_of_scope_stereo_is_rejected(smiles: str) -> None:
    assert not SMILESNNamer().name(smiles).success
