# IUPAC: P-65.1.2
# Layer: L2,L3,L5
"""Open-chain saturated monohydroxyalkanoic acids (carboxylic acid parent + hydroxy).

Carboxyl is principal characteristic group; non-carboxyl OH is hydroxy prefix
with locant from carboxyl numbering (COOH carbon = 1).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: α/β/ω-hydroxy monoacids
    ("OC(=O)C(O)C", "2-hydroxypropanoic acid", "2-羟基丙酸"),
    ("OC(=O)CCO", "3-hydroxypropanoic acid", "3-羟基丙酸"),
    ("OC(=O)CC(O)C", "3-hydroxybutanoic acid", "3-羟基丁酸"),
    ("O=C(O)CCCCCO", "6-hydroxyhexanoic acid", "6-羟基己酸"),
    ("OC(=O)C(O)CC", "2-hydroxybutanoic acid", "2-羟基丁酸"),
    # negative: plain acid, monoalcohol, diol must not become hydroxyacids
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCO", "ethanol", "乙醇"),
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hydroxyalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
