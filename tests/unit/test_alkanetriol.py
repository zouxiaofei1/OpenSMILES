# IUPAC: P-63.1.1
# Layer: L2,L4,L5
"""Unsubstituted open-chain saturated alkanetriols (exactly three OH).

Parent chain through all three OH carbons; name alkane-{a},{b},{c}-triol /
{烷}-{a},{b},{c}-三醇. Diols and monoalcohols must not break.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted open-chain triols
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("OCC(O)CCO", "butane-1,2,4-triol", "丁烷-1,2,4-三醇"),
    ("OCC(O)CCCO", "pentane-1,2,5-triol", "戊烷-1,2,5-三醇"),
    ("CC(O)C(O)CO", "butane-1,2,3-triol", "丁烷-1,2,3-三醇"),
    ("OCCC(O)CO", "butane-1,2,4-triol", "丁烷-1,2,4-三醇"),
    ("OCC(O)C(O)C", "butane-1,2,3-triol", "丁烷-1,2,3-三醇"),
    # negative: monoalcohol / diol / cycloalcohol must not break
    ("CCO", "ethanol", "乙醇"),
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
    ("OCCCO", "propane-1,3-diol", "丙烷-1,3-二醇"),
    ("C1CCC(O)CC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanetriol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
