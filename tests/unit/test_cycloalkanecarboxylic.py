# IUPAC: P-65.1.1
# Layer: L2,L4,L5
"""Monocyclic cycloalkanecarboxylic acids (unsub or mono ring sub).

Parent = saturated monocarbocycle with one exocyclic COOH on a ring carbon.
Allow 0–1 ring simple sub: monohalo (F/Cl/Br/I) or mono C1–C4 n-alkyl.
COOH attachment = locant 1; ring subs get lowest compatible locants.
IUPAC: cycloalkanecarboxylic acid / 环…烷甲酸.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted cycloalkanecarboxylic acids C3–C6
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),
    ("OC(=O)C1CCCC1", "cyclopentanecarboxylic acid", "环戊烷甲酸"),
    ("OC(=O)C1CCC1", "cyclobutanecarboxylic acid", "环丁烷甲酸"),
    ("OC(=O)C1CC1", "cyclopropanecarboxylic acid", "环丙烷甲酸"),
    # positive: mono C1–C4 n-alkyl on ring (COOH = 1)
    ("CC1CCCCC1C(=O)O", "2-methylcyclohexanecarboxylic acid", "2-甲基环己烷甲酸"),
    ("CCC1CCCCC1C(=O)O", "2-ethylcyclohexanecarboxylic acid", "2-乙基环己烷甲酸"),
    ("CC1CCCC1C(=O)O", "2-methylcyclopentanecarboxylic acid", "2-甲基环戊烷甲酸"),
    # positive: monohalo on ring
    ("ClC1CCCCC1C(=O)O", "2-chlorocyclohexanecarboxylic acid", "2-氯环己烷甲酸"),
    ("BrC1CCCCC1C(=O)O", "2-bromocyclohexanecarboxylic acid", "2-溴环己烷甲酸"),
    ("FC1CCCC1C(=O)O", "2-fluorocyclopentanecarboxylic acid", "2-氟环戊烷甲酸"),
    # negative: arene acid / open-chain acid / cycloalcohol must not change
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("OC(=O)CCCCCC", "heptanoic acid", "庚酸"),
    ("CCCCCCC(=O)O", "heptanoic acid", "庚酸"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkanecarboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cyclohexanecarboxylic_not_open_chain() -> None:
    """Ring COOH must not collapse into heptanoic acid."""
    r = SMILESNNamer().name("OC(=O)C1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexanecarboxylic acid"
    assert "heptanoic" not in en


def test_methyl_cycloacid_not_octanoic() -> None:
    """Monoalkyl ring acid must not collapse to open-chain octanoic acid."""
    r = SMILESNNamer().name("CC1CCCCC1C(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylcyclohexanecarboxylic acid"
    assert "octanoic" not in en
