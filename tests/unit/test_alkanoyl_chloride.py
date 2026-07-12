# IUPAC: P-65.5.1
# Layer: L1,L2,L3,L4,L5
"""Unsubstituted open-chain alkanoyl chlorides (mono acyl chlorides).

Exactly one –C(=O)Cl; parent chain through acyl carbon; no other main FG,
no ring/unsaturation this round. EN: C2 acetyl chloride; C>=3 …oyl chloride.
ZH: …酰氯 (乙酰氯, 丙酰氯, …).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: straight-chain mono alkanoyl chlorides
    ("ClC(=O)C", "acetyl chloride", "乙酰氯"),
    ("CC(=O)Cl", "acetyl chloride", "乙酰氯"),
    ("ClC(=O)CC", "propanoyl chloride", "丙酰氯"),
    ("ClC(=O)CCC", "butanoyl chloride", "丁酰氯"),
    ("ClC(=O)CCCC", "pentanoyl chloride", "戊酰氯"),
    ("ClC(=O)CCCCC", "hexanoyl chloride", "己酰氯"),
    # negative: acid / aldehyde / ketone / dichloroalkane must not become acyl chloride
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("ClCCCl", "1,2-dichloroethane", "1,2-二氯乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanoyl_chloride(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
