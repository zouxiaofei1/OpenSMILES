# IUPAC: P-29.3.1
# Layer: L3,L4,L5
"""Straight-chain monoalkyl (C1–C4) side-chain prefixes on alkane/alcohol parents."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: methyl on alcohol / alkane
    ("CC(C)CO", "2-methylpropan-1-ol", "2-甲基丙-1-醇"),
    ("CC(C)(C)O", "2-methylpropan-2-ol", "2-甲基丙-2-醇"),
    ("CCC(C)(C)O", "2-methylbutan-2-ol", None),
    ("CC(C)C", "2-methylpropane", "2-甲基丙烷"),
    ("CCC(C)C", "2-methylbutane", "2-甲基丁烷"),
    # negative: unsubstituted chain alcohols must stay correct
    ("CCCO", "propan-1-ol", "丙-1-醇"),
    ("CC(O)C", "propan-2-ol", "丙-2-醇"),
    ("CC(O)CC", "butan-2-ol", "丁-2-醇"),
    ("CCO", "ethanol", "乙醇"),
    ("CO", "methanol", "甲醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_monoalkyl_chain(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
