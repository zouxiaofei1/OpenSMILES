# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Simple mono free sulfonic acid / sulfonate anion parents (P-65.3).

Scope (first cut):
- open-chain alkanesulfonic acid R–SO2–OH, R = n-alkyl C1–C4
- free anion R–SO2–O⁻ → …sulfonate / …磺酸根 via parent anion flag
- alkali metal salt via L0 salt dissociation (sodium …sulfonate)
- unfused arenesulfonic acid Ar–SO2–OH, Ph with ≤2 Me/halo leaves

Out of scope: sulfonate esters (kind sulfonate), poly-acids, fused arenes,
vinyl/styrene leaves (no simple leaf), long chains >C4.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # unsubstituted alkanesulfonic acid
    ("CS(=O)(=O)O", "methanesulfonic acid", "甲磺酸"),
    ("CCS(=O)(=O)O", "ethanesulfonic acid", "乙磺酸"),
    ("CCCS(=O)(=O)O", "propanesulfonic acid", "丙磺酸"),
    ("CCCCS(=O)(=O)O", "butanesulfonic acid", "丁磺酸"),
    # free anion
    ("CS(=O)(=O)[O-]", "methanesulfonate", "甲磺酸根"),
    (
        "Cc1ccc(S(=O)(=O)[O-])cc1",
        "4-methylbenzenesulfonate",
        "4-甲基苯磺酸根",
    ),
    # alkali metal salt (L0 + anion + metal wrapper)
    ("CS(=O)(=O)[O-].[Na+]", "sodium methanesulfonate", "甲磺酸钠"),
    # unsubstituted / simple arenesulfonic acid
    ("c1ccc(S(=O)(=O)O)cc1", "benzenesulfonic acid", "苯磺酸"),
    (
        "Cc1ccc(S(=O)(=O)O)cc1",
        "4-methylbenzenesulfonic acid",
        "4-甲基苯磺酸",
    ),
    (
        "Clc1ccccc1S(=O)(=O)O",
        "2-chlorobenzenesulfonic acid",
        "2-氯苯磺酸",
    ),
    # negative near-miss: must keep existing correct names
    (
        "CC1=CC=C(C=C1)S(=O)(=O)OCCCC",
        "butyl 4-methylbenzenesulfonate",
        "对甲苯磺酸正丁酯",
    ),
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride", "甲磺酰氯"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


# Must not be named as free sulfonic acid / collapse to alkane
NEG_NOT_SULFONIC = [
    ("CC1=CC=C(C=C1)S(=O)(=O)OCCCC", "sulfonic acid"),
    ("c1ccc(S(=O)(=O)N)cc1", "sulfonic acid"),
    ("CS(=O)(=O)Cl", "sulfonic acid"),
    ("CC(=O)O", "sulfonic"),
    ("CS(=O)(=O)O", "methane"),  # must not collapse to methane alone
]


# First-cut aryl: long n-alkyl / vinyl leaves must not get wrong short names
NEG_ARYL_SCOPE = [
    (
        "CCCCCCCCCCCCc1ccc(S(=O)(=O)O)cc1",
        "4-butylbenzenesulfonic acid",
    ),
    (
        "[Na+].C=CC1=CC=C(C=C1)S(=O)(=O)[O-]",
        "sodium 4-ethylbenzenesulfonate",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_sulfonic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONIC)
def test_not_wrong_sulfonic(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    if forbidden == "methane":
        assert en != "methane"
        assert "sulfonic" in en or "sulfonate" in en
        return
    if forbidden == "sulfonic":
        # carboxylic must not pick up sulfonic wording
        assert "sulfonic" not in en
        assert "磺" not in normalize_zh(r.zh)
        return
    assert forbidden not in en
    assert "磺酸" not in normalize_zh(r.zh) or "酯" in normalize_zh(r.zh) or "酰胺" in normalize_zh(r.zh) or "酰氯" in normalize_zh(r.zh)


@pytest.mark.parametrize("smiles,wrong", NEG_ARYL_SCOPE)
def test_aryl_first_cut_scope(smiles: str, wrong: str) -> None:
    """Long n-alkyl / vinyl Ph must not produce truncated Me/Et/Bu wrong names."""
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en) if r.success else ""
    assert normalize_en(wrong) not in en
    assert "butylbenzene" not in en
    assert "ethylbenzene" not in en


def test_l1_detects_free_acid_not_ester() -> None:
    mol = preprocess("CS(=O)(=O)O")
    assert mol is not None
    info = analyze(mol)
    assert info.get("has_sulfonic_acid")
    assert len(info.get("sulfonic_acids") or []) == 1
    assert not info.get("has_sulfonate")


def test_l1_detects_anion() -> None:
    mol = preprocess("CS(=O)(=O)[O-]")
    assert mol is not None
    info = analyze(mol)
    assert info.get("has_sulfonic_acid")
    e = (info.get("sulfonic_acids") or [None])[0]
    assert e is not None and e.get("anion") is True


def test_l1_ester_not_free_acid() -> None:
    mol = preprocess("CC1=CC=C(C=C1)S(=O)(=O)OCCCC")
    assert mol is not None
    info = analyze(mol)
    assert info.get("has_sulfonate")
    assert not info.get("has_sulfonic_acid")


def test_methanesulfonic_not_methane() -> None:
    r = SMILESNNamer().name("CS(=O)(=O)O")
    assert r.success
    assert normalize_en(r.en) == "methanesulfonic acid"
    assert normalize_en(r.en) != "methane"
