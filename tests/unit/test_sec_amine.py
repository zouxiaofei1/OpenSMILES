# IUPAC: P-62.2.2.1
# Layer: L1,L2,L3,L5
"""Open-chain secondary monoamines: N-alkylalkanamines (symmetric preferred).

Exactly one secondary amine N (two C neighbors, one H); both N-substituents
unsubstituted linear alkyl C1–C4. Parent = longer alkyl as alkanamine;
shorter as N-alkyl prefix. English: N-ethylethanamine (no space after N-).
Chinese: N-乙基乙胺.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: symmetric and simple unsymmetric secondary amines
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("CCNC", "N-methylethanamine", "N-甲基乙胺"),
    ("CNC", "N-methylmethanamine", "N-甲基甲胺"),
    ("CCCNCCC", "N-propylpropan-1-amine", "N-丙基丙-1-胺"),
    ("CCNCCC", "N-ethylpropan-1-amine", "N-乙基丙-1-胺"),
    ("CCCNC", "N-methylpropan-1-amine", "N-甲基丙-1-胺"),
    # negative: primary monoamine / diamine / tertiary must not break
    ("CCN", "ethanamine", "乙胺"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("CCN(CC)CC", "N,N-diethylethanamine", "N,N-二乙基乙胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sec_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_primary_amine_form() -> None:
    """Secondary CCNCC must not collapse to ethane or ethanamine."""
    r = SMILESNNamer().name("CCNCC")
    assert r.success
    assert normalize_en(r.en) == "n-ethylethanamine"
    assert "ethane" != normalize_en(r.en)
