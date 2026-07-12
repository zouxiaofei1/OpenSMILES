# IUPAC: P-62.2.1
# Layer: L2,L4,L5
"""Unsubstituted open-chain saturated alkanediamines (exactly two primary amines).

Parent chain through both amine carbons; name alkane-{a},{b}-diamine /
{烷}-{a},{b}-二胺. Side chains use existing L3 alkyl prefixes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted open-chain diamines
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("NCCCN", "propane-1,3-diamine", "丙烷-1,3-二胺"),
    ("NCCCCN", "butane-1,4-diamine", "丁烷-1,4-二胺"),
    ("NCCCCCN", "pentane-1,5-diamine", "戊烷-1,5-二胺"),
    ("NC(C)CCN", "butane-1,3-diamine", "丁烷-1,3-二胺"),
    ("NCC(C)CN", "2-methylpropane-1,3-diamine", "2-甲基丙烷-1,3-二胺"),
    # negative: monoamines / cycloamine must not break
    ("CCN", "ethanamine", "乙胺"),
    ("CC(N)C", "propan-2-amine", "丙-2-胺"),
    ("C1CCC(N)CC1", "cyclohexanamine", "环己胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanediamine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
