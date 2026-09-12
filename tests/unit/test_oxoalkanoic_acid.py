# IUPAC: P-65.1.2
# Layer: L3
"""Open-chain saturated mono-oxoalkanoic acids (carboxylic acid parent + oxo).

Carboxyl is principal characteristic group; ketone carbonyl is oxo/氧代 prefix
with locant from carboxyl numbering (COOH carbon = 1).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: α/β/γ-oxo monoacids
    ("OC(=O)C(=O)C", "2-oxopropanoic acid", "2-氧代丙酸"),
    ("OC(=O)CCC(=O)C", "4-oxopentanoic acid", "4-氧代戊酸"),
    ("CC(=O)CC(=O)O", "3-oxobutanoic acid", "3-氧代丁酸"),
    ("CCCC(=O)C(=O)O", "2-oxopentanoic acid", None),
    # negative: plain acid, ketone, dione, amino acid must not become oxoacids
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
    ("NCC(=O)O", "2-aminoacetic acid", "2-氨基乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_oxoalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
