# IUPAC: P-63.1.5
# Layer: L1,L2,L4,L5
"""Unsubstituted open-chain monovalent thiols (alkanethiols).

Exactly one –SH (S with one C neighbor and ≥1 H); no other principal FG;
saturated acyclic. Parent chain through SH-attached carbon; locants like alcohol.
C1–C2 omit locant (methanethiol/ethanethiol); C≥3: alkane-n-thiol / 首字-n-硫醇.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # C2 烯硫醇：烯 1-2 + SH 1 无歧义，位次省略融合 ethenethiol (P-14.3.4)
    ("C=CS", "ethenethiol", "乙烯硫醇"),
    # alkanethiol 正例（P-63.1.5）：C1–C2 omit locant，C≥3 用 alkane-n-thiol
    ("CS", "methanethiol", "甲硫醇"),
    ("CCS", "ethanethiol", "乙硫醇"),
    ("CCCS", "propane-1-thiol", "丙-1-硫醇"),
    ("CCCCCS", "pentane-1-thiol", "戊-1-硫醇"),
    ("CC(C)S", "propane-2-thiol", "丙-2-硫醇"),
    ("CCCC(CC)S", "hexane-3-thiol", "己-3-硫醇"),
    # 二硫醇（数量后缀生成式，对齐 amine）
    ("SCCS", "ethane-1,2-dithiol", "乙烷-1,2-二硫醇"),
    ("SCCCCS", "butane-1,4-dithiol", "丁烷-1,4-二硫醇"),
    # negative: alcohol / amine / thioether (not thiol)
    ("CCO", "ethanol", "乙醇"),
    ("CCN", "ethanamine", "乙胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanethiol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
