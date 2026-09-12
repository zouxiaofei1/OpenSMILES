# IUPAC: P-66.6.1 / P-65.6 / P-44
# Layer: L2
"""Reject formaldehyde / methyl formate false parents for ring-only FG carbons.

Aldehyde or ester carbonyls that only attach to ring/aryl carbons must not
collapse to open-chain C1 formaldehyde or methyl formate. True H2C=O and
HCO2Me remain.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

POS_NO_FORMALDEHYDE = [
    "O=Cc1ccc(F)c(F)c1COc1ccccc1",
    "O=Cc1ccncc1",
]

POS_NO_METHYL_FORMATE = [
    "COC(=O)c1cccs1",  # methyl thiophene-2-carboxylate-like
]

NEG = [
    ("C=O", "formaldehyde", "甲醛"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("COC=O", "methyl formate", "甲酸甲酯"),
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
]


@pytest.mark.parametrize("smiles", POS_NO_FORMALDEHYDE)
def test_ring_aldehyde_not_formaldehyde(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "formaldehyde"
    assert "formaldehyde" not in en


@pytest.mark.parametrize("smiles", POS_NO_METHYL_FORMATE)
def test_ring_ester_not_methyl_formate(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methyl formate"


@pytest.mark.parametrize("smiles,en,zh", NEG)
def test_formyl_formate_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
