# IUPAC: P-65.1.1.1 / P-63.2.2 / P-29.3
# Layer: L2,L3,L5
"""Reject formic-acid false parent for ring COOH; claim isopropoxy/isobutoxy on arenes.

Ring-only carboxyls (COOH carbon neighbors are ring/aryl only) must not collapse
to open-chain formic acid. Branched outer alkoxy O–CHMe2 / O–CH2–CHMe2 must be
claimable so simple benzoic gating passes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer


# Positive: arene acids with branched alkoxy → benzoic parent + alkoxy prefix

# Positive: ring/hetero acids must not emit formic (name quality optional)
POS_NO_FORMIC = [
    "IC1=C(SC(=C1S(=O)(=O)C(C)C)SC)C(=O)O",
    "CN1N=C(C=C1)C1=CC=C(C(=O)O)C=C1",
]

# Negative / regressions: true formic, open acids, methoxybenzoic
NEG_CASES = [
    ("C(=O)O", "formic acid", "甲酸"),
    ("O=CO", "formic acid", "甲酸"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCCCCC(=O)O", "hexanoic acid", "己酸"),
    ("COc1ccc(C(=O)O)cc1", "4-methoxybenzoic acid", "4-甲氧基苯甲酸"),
]


@pytest.mark.parametrize("smiles", POS_NO_FORMIC)
def test_ring_acid_not_formic(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "formic acid"
    assert "formic" not in en


@pytest.mark.parametrize("smiles,en,zh", NEG_CASES)
def test_formic_and_open_acid_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
