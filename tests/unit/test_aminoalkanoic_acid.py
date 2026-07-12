# IUPAC: P-65.1.2
# Layer: L3,L5
"""Open-chain saturated monoaminoalkanoic acids (carboxylic acid parent + amino).

Carboxyl is principal characteristic group; primary amino is amino prefix
with locant from carboxyl numbering (COOH carbon = 1).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: α/β/ω-amino monoacids
    ("NCC(=O)O", "2-aminoacetic acid", "2-氨基乙酸"),
    ("OC(=O)C(N)C", "2-aminopropanoic acid", "2-氨基丙酸"),
    ("NCCC(=O)O", "3-aminopropanoic acid", "3-氨基丙酸"),
    ("CCCCC(N)C(=O)O", "2-aminohexanoic acid", None),
    ("NCCCC(=O)O", "4-aminobutanoic acid", "4-氨基丁酸"),
    # negative: plain acid, monoamine, hydroxyacid must not break
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCN", "ethanamine", "乙胺"),
    ("OC(=O)C(O)C", "2-hydroxypropanoic acid", "2-羟基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_aminoalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
