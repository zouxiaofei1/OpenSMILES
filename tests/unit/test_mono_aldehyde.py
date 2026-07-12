# IUPAC: P-66.6.1
# Layer: L1,L2,L4,L5
"""Simple acyclic monoaldehydes (alkanals): retained C1/C2 + systematic C3+.

Aldehyde carbonyl carbon has exactly one carbon neighbor and is not carboxyl.
Suffix -al / 醛 has no locant (aldehyde carbon is always position 1).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: retained + straight-chain + branched monoaldehydes
    ("C=O", "formaldehyde", "甲醛"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCC=O", "propanal", "丙醛"),
    ("CCCCC=O", "pentanal", "戊醛"),
    ("CCCCCC=O", "hexanal", "己醛"),
    ("CC(C)C=O", "2-methylpropanal", "2-甲基丙醛"),
    # negative: alkane / ketone / acid / alcohol must not become aldehydes
    ("CCC", "propane", "丙烷"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_aldehyde(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
