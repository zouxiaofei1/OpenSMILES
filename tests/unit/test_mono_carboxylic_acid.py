# IUPAC: P-65.1.1
# Layer: L1,L2,L4,L5
"""Monocarboxylic acids: systematic oic acid + retained formic/acetic.

Carboxyl carbon is part of the parent chain (locant 1); suffix has no locant.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: retained + straight-chain + branched monoacids
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCC(=O)O", "propanoic acid", "丙酸"),
    ("CCCC(=O)O", "butanoic acid", "丁酸"),
    ("CCCCC(=O)O", "pentanoic acid", "戊酸"),
    ("CCCCCC(=O)O", "hexanoic acid", "己酸"),
    ("CC(C)C(=O)O", "2-methylpropanoic acid", "2-甲基丙酸"),
    ("CC(C)CCC(=O)O", "4-methylpentanoic acid", "4-甲基戊酸"),
    ("C(=O)O", "formic acid", "甲酸"),
    # negative: alcohols and haloalkanes must not become acids
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
    ("CCCl", "chloroethane", "氯乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_carboxylic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
