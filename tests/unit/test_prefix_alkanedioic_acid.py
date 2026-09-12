# IUPAC: P-65.1.2
# Layer: L2,L3,L4,L5
"""Open-chain saturated substituted alkanedioic acids / dioates.

Parent remains alkanedioic acid / …dioate (or oxalic); chain hydroxy /
amino / oxo are prefixes (incl. multi-prefix); numbering from either
carboxyl carbon by lowest set. Anion → …dioate / …二酸根.
"""
from __future__ import annotations

import pytest

from rdkit import Chem

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: hydroxy / amino / oxo prefixes on saturated diacids
    ("O=C(O)CC(O)CC(=O)O", "3-hydroxypentanedioic acid", "3-羟基戊二酸"),
    ("NC(CC(=O)O)CC(=O)O", "3-aminopentanedioic acid", "3-氨基戊二酸"),
    ("O=C(O)C(N)CC(=O)O", "2-aminobutanedioic acid", "2-氨基丁二酸"),
    ("O=C(O)C(O)C(=O)O", "2-hydroxypropanedioic acid", "2-羟基丙二酸"),
    ("O=C([O-])CCCC(=O)C(=O)[O-]", "2-oxohexanedioate", "2-氧代己二酸根"),
    ("O=C(O)C(=O)C(O)C(=O)O", "2-hydroxy-3-oxobutanedioic acid", "2-羟基-3-氧代丁二酸"),
    ("O=C(O)C(O)C(O)C(O)C(=O)O", "2,3,4-trihydroxypentanedioic acid", "2,3,4-三羟基戊二酸"),
    ("O=C(O)C(C)(C)C(O)C(=O)O", "3-hydroxy-2,2-dimethylbutanedioic acid", "3-羟基-2,2-二甲基丁二酸"),
    # negative: unsubstituted diacid, mono hydroxyacid, alkyl-diacid, monoacid
    ("O=C(O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("O=C(O)C(O)C", "2-hydroxypropanoic acid", "2-羟基丙酸"),
    ("O=C(O)C(C)C(=O)O", "2-methylpropanedioic acid", "2-甲基丙二酸"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_prefix_alkanedioic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
