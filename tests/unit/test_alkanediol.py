# IUPAC: P-63.1.1
# Layer: L2,L4,L5
"""Unsubstituted open-chain saturated alkanediols (exactly two OH).

Parent chain through both OH carbons; name alkane-{a},{b}-diol /
{烷}-{a},{b}-二醇. Side chains use existing L3 alkyl prefixes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted open-chain diols
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
    ("OCCCO", "propane-1,3-diol", "丙烷-1,3-二醇"),
    ("CC(O)C(C)O", "butane-2,3-diol", "丁烷-2,3-二醇"),
    ("OCC(C)O", "propane-1,2-diol", "丙烷-1,2-二醇"),
    ("CC(C)C(O)CO", "3-methylbutane-1,2-diol", "3-甲基丁烷-1,2-二醇"),
    ("OCCCCO", "butane-1,4-diol", "丁烷-1,4-二醇"),
    # negative: monoalcohols / cycloalcohol must not break
    ("CCO", "ethanol", "乙醇"),
    ("CC(O)C", "propan-2-ol", "丙-2-醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanediol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
