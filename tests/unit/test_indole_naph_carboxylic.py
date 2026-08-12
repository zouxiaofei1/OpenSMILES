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
