# IUPAC: P-65.1.1 / P-72.2.2.1
# Layer: L1,L2,L5
"""Carboxylate anions: C(=O)[O-] → alkanoate / …酸根 (not aldehyde).

IUPAC P-65.1.1 parent acid chain numbering; P-72.2.2.1 anion suffix -ate / 酸根.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: bare carboxylate anions
    ("CCCCCCCCCCCC(=O)[O-]", "dodecanoate", "十二酸根"),
    ("CCCCCCCCCCCCCC(=O)[O-]", "tetradecanoate", "十四酸根"),
    ("CC(C)CC(=O)[O-]", "3-methylbutanoate", "3-甲基丁酸根"),
    ("CC(=O)[O-]", "acetate", "乙酸根"),
    ("CCCCCCCCCCC(O)C(=O)[O-]", "2-hydroxydodecanoate", "2-羟基十二酸根"),
    ("NCCCC(=O)[O-]", "4-aminobutanoate", "4-氨基丁酸根"),
    ("C(=O)[O-]", "formate", "甲酸根"),
    # negative: neutral acid, true aldehyde, long alkane must not regress
    ("CCCCCCCCCCCC(=O)O", "dodecanoic acid", "十二酸"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCC", "hexacosane", "二十六烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_carboxylate_anion(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
