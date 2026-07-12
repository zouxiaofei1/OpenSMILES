# IUPAC: P-65.1.1
# Layer: L2,L5
"""Unsubstituted open-chain saturated alkanedioic acids (exactly two COOH).

Parent chain through both carboxyl carbons; name alkanedioic acid /
{烷首}二酸. C2 retained: oxalic acid / 草酸. No locants for terminal diacids.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted open-chain diacids
    ("OC(=O)C(=O)O", "oxalic acid", "草酸"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("OC(=O)CCC(=O)O", "butanedioic acid", "丁二酸"),
    ("OC(=O)CCCC(=O)O", "pentanedioic acid", "戊二酸"),
    ("O=C(O)CCCCCCCC(=O)O", "nonanedioic acid", "壬二酸"),
    # negative: monoacids must not become diacids
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCC(=O)O", "propanoic acid", "丙酸"),
    ("CCCC(=O)O", "butanoic acid", "丁酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanedioic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
