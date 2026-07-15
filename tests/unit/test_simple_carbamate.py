# IUPAC: P-65
# Layer: L1,L2,L5
"""Simple alkyl N-…carbamate (incl. Boc / tert-butyl carbamate).

Carbamate R2N–C(=O)–OR is principal FG (not ester). English follows gold:
  methyl N-methylcarbamate; tert-butyl (3-methoxyphenyl)carbamate.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: simple alkyl N-alkyl carbamates
    ("CNC(=O)OC", "methyl N-methylcarbamate", "N-甲基氨基甲酸甲酯"),
    ("CCNC(=O)OCC", "ethyl N-ethylcarbamate", "N-乙基氨基甲酸乙酯"),
    # positive: Boc / tert-butyl (gold style without forced N- for aryl)
    (
        "CC(C)(C)OC(=O)Nc1ccccc1",
        "tert-butyl N-phenylcarbamate",
        "N-苯基氨基甲酸叔丁酯",
    ),
    (
        "CC(C)(C)OC(=O)NCc1ccccc1",
        "tert-butyl N-benzylcarbamate",
        "N-苄基氨基甲酸叔丁酯",
    ),
    (
        "COC=1C=C(C=CC1)NC(OC(C)(C)C)=O",
        "tert-butyl (3-methoxyphenyl)carbamate",
        None,
    ),
    (
        "CN(C(OC(C)(C)C)=O)C1=CC(=CC=C1)[N+](=O)[O-]",
        "tert-butyl methyl(3-nitrophenyl)carbamate",
        None,
    ),
    # negative: true ester / amide must not become carbamate
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_carbamate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
