# IUPAC: P-68.1
# Layer: L1,L2,L4,L5
"""Simple arylboronic acids (Ar–B(OH)2): phenylboronic acid parent (P-68.1).

Unsubstituted phenylboronic acid and ring-substituted halo/methyl/CF3
variants with B-attachment as locant 1. Must not steal carboxylic acids,
phenol, benzene, or sulfoxides.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("B(O)(O)c1ccccc1", "phenylboronic acid", "苯硼酸"),
    ("c1ccccc1B(O)O", "phenylboronic acid", "苯硼酸"),
    # positive: mono-sub methyl
    ("Cc1ccc(B(O)O)cc1", "4-methylphenylboronic acid", "4-甲基苯基硼酸"),
    # positive: dual gold (Br + CF3)
    (
        "BrC=1C=C(C=C(C1)C(F)(F)F)B(O)O",
        "3-bromo-5-(trifluoromethyl)phenylboronic acid",
        "3-溴-5-三氟甲基苯基硼酸",
    ),
    # negative: must not mis-claim
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccccc1O", "phenol", "苯酚"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CS(C)=O", "dimethyl sulfoxide", "二甲基亚砜"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_boronic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
