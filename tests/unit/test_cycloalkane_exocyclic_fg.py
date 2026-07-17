# IUPAC: P-66.6.1 / P-66.5.1 / P-66.1.1 / P-65.6 / P-65.5
# Layer: L2,L4,L5
"""Monocyclic cycloalkane with one exocyclic carbonyl-class FG (P-65/P-66).

Parent = unfused saturated monocarbocycle C3–C10 with exactly one ring-attached
exocyclic FG of type: carbaldehyde / carbonitrile / carboxamide /
carboxylate (alkyl ester) / carbonyl halide.
Ring may carry 0–1 simple sub: monohalo or mono C1–C4 n-alkyl; FG attach = 1.
IUPAC PIN stems: cyclohexanecarbaldehyde / carbonitrile / carboxamide /
carboxylate / carbonyl chloride (aligned with cycloalkanecarboxylic acid).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted carbaldehyde / carbonitrile / carboxamide
    ("O=CC1CCCCC1", "cyclohexanecarbaldehyde", "环己烷甲醛"),
    ("O=CC1CCCC1", "cyclopentanecarbaldehyde", "环戊烷甲醛"),
    ("N#CC1CCCCC1", "cyclohexanecarbonitrile", "环己烷甲腈"),
    ("N#CC1CCCC1", "cyclopentanecarbonitrile", "环戊烷甲腈"),
    ("NC(=O)C1CCCCC1", "cyclohexanecarboxamide", "环己烷甲酰胺"),
    ("NC(=O)C1CCCC1", "cyclopentanecarboxamide", "环戊烷甲酰胺"),
    # positive: simple alkyl ester / acyl chloride (same ring-FG topology)
    ("COC(=O)C1CCCCC1", "methyl cyclohexanecarboxylate", "环己烷甲酸甲酯"),
    ("ClC(=O)C1CCCCC1", "cyclohexanecarbonyl chloride", "环己烷甲酰氯"),
    # positive: mono ring sub (FG attach = 1)
    ("CC1CCCCC1C=O", "2-methylcyclohexanecarbaldehyde", "2-甲基环己烷甲醛"),
    ("ClC1CCCCC1C=O", "2-chlorocyclohexanecarbaldehyde", "2-氯环己烷甲醛"),
    # negative: open-chain formic-style must stay; arene / cyclo acid unchanged
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC#N", "acetonitrile", "乙腈"),
    ("CC(N)=O", "acetamide", "乙酰胺"),
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkane_exocyclic_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cyclohexanecarbaldehyde_not_ethane() -> None:
    """Ring CHO must not collapse to formic-parent ethane."""
    r = SMILESNNamer().name("O=CC1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexanecarbaldehyde"
    assert "ethane" not in en


def test_cyclohexanecarbonitrile_not_formonitrile() -> None:
    """Ring CN must not use formonitrile + cyclohexyl sub."""
    r = SMILESNNamer().name("N#CC1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexanecarbonitrile"
    assert "formonitrile" not in en
    assert "cyclohexyl" not in en
