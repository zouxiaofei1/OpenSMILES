# IUPAC: P-63.1.1
# Layer: L2,L4,L5
"""Open-chain saturated polyols beyond tetraol (multiplicity-generic).

The count suffix (diol/triol/tetraol/pentaol/…) is generated as MULT + base
suffix, so any multiplicity is supported without hard-coded variant entries.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("OCC(O)C(O)C(O)CO", "pentane-1,2,3,4,5-pentaol", "戊烷-1,2,3,4,5-五醇"),
    ("OCC(O)C(O)C(O)C(O)CO", "hexane-1,2,3,4,5,6-hexaol", "己烷-1,2,3,4,5,6-六醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyol_extension(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
