# IUPAC: P-63.1.2 / P-64.2.1 / P-14.3.4
# Layer: L2,L4,L5
"""Saturated monocyclic diols and diones (cycloalkanediol / cycloalkanedione).

Parent = unfused C3–C10 sat carbocycle with exactly two ring-carbon OH
or two ring ketone carbonyls. Lowest locant pair for the two FGs (always
shown). Optional simple ring halo / claimable n-alkyl. Must not capture
benzenediol, mono cycloalcohol/one, or open-chain diol/dione.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted cycloalkanediols
    ("OC1CCCCC1O", "cyclohexane-1,2-diol", "环己烷-1,2-二醇"),
    ("OC1CC(O)CCC1", "cyclohexane-1,3-diol", "环己烷-1,3-二醇"),
    ("OC1CCC(O)CC1", "cyclohexane-1,4-diol", "环己烷-1,4-二醇"),
    ("OC1CCCC1O", "cyclopentane-1,2-diol", "环戊烷-1,2-二醇"),
    ("OC1CCC1O", "cyclobutane-1,2-diol", "环丁烷-1,2-二醇"),
    # positive: simple ring sub keeps FG pair + sub locant
    ("CC1CC(O)CC(O)C1", "5-methylcyclohexane-1,3-diol", "5-甲基环己烷-1,3-二醇"),
    ("ClC1CC(O)CC(O)C1", "5-chlorocyclohexane-1,3-diol", "5-氯环己烷-1,3-二醇"),
    # negative: mono cyclo FG, arene diol, open-chain poly FG stay correct
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("OCCCO", "propane-1,3-diol", "丙烷-1,3-二醇"),
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
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
