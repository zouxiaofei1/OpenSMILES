# IUPAC: P-63.1.2
# Layer: L2,L3,L4,L5
"""Unsubstituted benzene-a,b-diol (exactly two phenolic OH on benzene).

Systematic names only: benzene-1,2-diol / benzene-1,3-diol / benzene-1,4-diol.
Chinese uses 二酚 (gold), not 二醇. No retained catechol/resorcinol/hydroquinone.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
# zh=None skips Chinese assert (chebi rows without chinese gold).
CASES = [
    # negative: mono phenol, chain diol/triol, bare benzene, chain alcohol
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("c1ccccc1", "benzene", "苯"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzenediol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
