# IUPAC: P-62.2.1
# Layer: L2,L4,L5
"""Open-chain saturated polyamines beyond tetraamine (multiplicity-generic).

Count suffix (diamine/triamine/tetraamine/pentaamine/…) is generated as
MULT + base suffix; no hard-coded variant entries.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("NCC(N)C(N)C(N)CN", "pentane-1,2,3,4,5-pentaamine", "戊烷-1,2,3,4,5-五胺"),
    ("NCC(N)C(N)C(N)C(N)CN", "hexane-1,2,3,4,5,6-hexaamine", "己烷-1,2,3,4,5,6-六胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyamine_extension(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
