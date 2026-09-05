# IUPAC: P-62.2.1
# Layer: L2,L3,L4,L5
"""Unsubstituted open-chain saturated polyamines (triamine, tetraamine).

Parent chain through all amine carbons; name alkane-{a},{b},{c}-triamine /
{烷}-{a},{b},{c}-三胺 (triamine) or alkane-{a},{b},{c},{d}-tetraamine /
{烷}-{a},{b},{c},{d}-四胺 (tetraamine).

Diamines and monoamines must not break.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # triamine: unsubstituted open-chain
    ("NCC(N)CN", "propane-1,2,3-triamine", "丙烷-1,2,3-三胺"),
    ("NCC(N)CCN", "butane-1,2,4-triamine", "丁烷-1,2,4-三胺"),
    ("NCCC(N)CCN", "pentane-1,3,5-triamine", "戊烷-1,3,5-三胺"),
    ("NCC(N)(C)CN", "2-methylpropane-1,2,3-triamine", "2-甲基丙烷-1,2,3-三胺"),
    ("NCC(N)C(N)C", "butane-1,2,3-triamine", "丁烷-1,2,3-三胺"),
    # multi-substituted triamine edge cases
    ("NCC(N)(CC)CN", "2-ethylpropane-1,2,3-triamine", "2-乙基丙烷-1,2,3-三胺"),
    ("NCC(C)C(N)CN", "3-methylbutane-1,2,4-triamine", "3-甲基丁烷-1,2,4-三胺"),
    # tetraamine: unsubstituted open-chain
    # negative: diamine / monoamine / triol must not break
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("NCCCN", "propane-1,3-diamine", "丙烷-1,3-二胺"),
    ("CCN", "ethanamine", "乙胺"),
    ("CC(N)C", "propan-2-amine", "丙-2-胺"),
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
