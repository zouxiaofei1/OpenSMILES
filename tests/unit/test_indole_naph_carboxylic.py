# IUPAC: P-65.1.1 / P-22.2.1 / P-25
# Layer: L2,L4,L5
"""Retained fused-ring carboxylic parents: indolecarboxylic / naphthalenecarboxylic.

IUPAC P-65.1.1 / P-22.2.1 / P-25: mono COOH on retained indole or naphthalene
core yields 1H-indole-n-carboxylic acid / naphthalene-n-carboxylic acid
(Chinese: 1H-吲哚-n-甲酸 / 萘-n-甲酸). Must not collapse to open-chain formic acid.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted indole-3-carboxylic
    (
        "O=C(O)c1c[nH]c2ccccc12",
        "1H-indole-3-carboxylic acid",
        "1H-吲哚-3-甲酸",
    ),
    # positive: 3,4-dimethyl-1H-indole-2-carboxylic
    (
        "Cc1cccc2[nH]c(C(=O)O)c(C)c12",
        "3,4-dimethyl-1H-indole-2-carboxylic acid",
        "3,4-二甲基-1H-吲哚-2-甲酸",
    ),
    # positive: 4,5-dimethoxy-1H-indole-2-carboxylic (zh gold)
    (
        "COC1=C2C=C(NC2=CC=C1OC)C(=O)O",
        "4,5-dimethoxy-1H-indole-2-carboxylic acid",
        "4,5-二甲氧基-1H-吲哚-2-甲酸",
    ),
    # positive: 5-methylnaphthalene-1-carboxylic
    (
        "Cc1cccc2c(C(=O)O)cccc12",
        "5-methylnaphthalene-1-carboxylic acid",
        "5-甲基萘-1-甲酸",
    ),
    # positive: unsubstituted naphthalene-2-carboxylic
    (
        "O=C(O)c1ccc2ccccc2c1",
        "naphthalene-2-carboxylic acid",
        "萘-2-甲酸",
    ),
    # positive: unsubstituted naphthalene-1-carboxylic
    (
        "O=C(O)c1cccc2ccccc12",
        "naphthalene-1-carboxylic acid",
        "萘-1-甲酸",
    ),
    # regression: benzoic still retained
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    # negative: open-chain acids must not become fused carboxylic parents
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    ("O=C(O)[C@@H](O)CO", "(2S)-2,3-dihydroxypropanoic acid", "(2S)-2,3-二羟基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indole_naph_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_indole_3_carboxylic_not_formic() -> None:
    """COOH on indole must not collapse into formic acid."""
    r = SMILESNNamer().name("O=C(O)c1c[nH]c2ccccc12")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1h-indole-3-carboxylic acid"
    assert "formic" not in en
