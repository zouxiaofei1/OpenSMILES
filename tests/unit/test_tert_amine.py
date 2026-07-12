# IUPAC: P-62.2.2.1
# Layer: L1,L2,L3,L5
"""Open-chain tertiary monoamines: N,N-dialkylalkanamines.

Exactly one tertiary amine N (three C neighbors, zero H); all three N-substituents
unsubstituted linear alkyl C1–C4. Parent = longest alkyl as alkanamine; the other
two as N-substituents. Same pair → N,N-di{alkyl}; different → N-{a}-N-{b}
alphabetical by alkyl English name. Chinese: N,N-二乙基乙胺.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: symmetric / near-symmetric tertiary amines
    ("CCN(CC)CC", "N,N-diethylethanamine", "N,N-二乙基乙胺"),
    ("CN(C)C", "N,N-dimethylmethanamine", "N,N-二甲基甲胺"),
    ("CCN(C)CC", "N-ethyl-N-methylethanamine", "N-乙基-N-甲基乙胺"),
    ("CCCN(CCC)CCC", "N,N-dipropylpropan-1-amine", "N,N-二丙基丙-1-胺"),
    ("CCN(C)C", "N,N-dimethylethanamine", "N,N-二甲基乙胺"),
    ("CCCN(C)C", "N,N-dimethylpropan-1-amine", "N,N-二甲基丙-1-胺"),
    # negative: primary / secondary / diamine must not break
    ("CCN", "ethanamine", "乙胺"),
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_tert_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_alkane_or_primary() -> None:
    """Tertiary CCN(CC)CC must not collapse to ethane or ethanamine."""
    r = SMILESNNamer().name("CCN(CC)CC")
    assert r.success
    assert normalize_en(r.en) == "n,n-diethylethanamine"
    assert "ethane" != normalize_en(r.en)
