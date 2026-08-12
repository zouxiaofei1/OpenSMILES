# IUPAC: P-62.2.2.1
# Layer: L1,L2,L3,L5
"""Open-chain tertiary monoamine scope (P-62.2.2.1) — negative guard only.

Positive tertiary-amine cases are not covered in this file. The retained
cases assert that a primary monoamine (CCN) and a diamine (NCCN) must not be
named as N,N-dialkylalkanamines.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: primary / secondary / diamine must not break
    ("CCN", "ethanamine", "乙胺"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_tert_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
