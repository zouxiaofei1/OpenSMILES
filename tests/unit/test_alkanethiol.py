# IUPAC: P-63.1.5
# Layer: L1,L2,L4,L5
"""Unsubstituted open-chain monovalent thiols (alkanethiols).

Exactly one –SH (S with one C neighbor and ≥1 H); no other principal FG;
saturated acyclic. Parent chain through SH-attached carbon; locants like alcohol.
C1–C2 omit locant (methanethiol/ethanethiol); C≥3: alkane-n-thiol / 首字-n-硫醇.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: simple alkanethiols
    ("CS", "methanethiol", "甲硫醇"),
    ("CCS", "ethanethiol", "乙硫醇"),
    ("CCCS", "propane-1-thiol", "丙-1-硫醇"),
    ("CC(C)S", "propane-2-thiol", "丙-2-硫醇"),
    ("CCCCS", "butane-1-thiol", "丁-1-硫醇"),
    # negative: alcohol / amine / thioether (keep current or non-thiol)
    ("CCO", "ethanol", "乙醇"),
    ("CCN", "ethanamine", "乙胺"),
    ("CCSC", "ethane", "乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanethiol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
