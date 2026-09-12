# IUPAC: P-66.1.1 / P-29 / P-22.2.1
# Layer: L2,L3,L5
"""N-heteroaryl benzamide via recursive substituent path."""
from __future__ import annotations

import pytest
from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer


REGRESS = [
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", REGRESS)
def test_n_heteroaryl_benzamide_regress(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
