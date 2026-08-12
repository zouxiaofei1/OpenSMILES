# IUPAC: P-25 / P-64
# Layer: L2,L4,L5
"""9,10-Anthraquinone retained parent (linear anthracene + meso dione).

PIN-style general name: 9,10-anthraquinone / 蒽醌 with ring methyl prefixes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # open-chain alkanedione not captured as anthraquinone
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
    # benzene untouched
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_anthraquinone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
