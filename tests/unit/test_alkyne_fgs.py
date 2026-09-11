# IUPAC: P-31.1 / P-63.1.1 / P-66.6.1 / P-66.5.1 / P-66.1.1 / P-14.3.4
# Layer: L4,L5
"""Open-chain mono alkyne FG parents: alkynol / alkynal / alkynenitrile / alkynamide.

Mono principal FG + mono C≡C (no C=C); FG carbon and triple-bond ends acyclic.
Kinds stay alcohol/aldehyde/nitrile/amide; unsaturation via triple_bond + yne_locant.

P-14.3.4：三碳炔的融合式后缀（-ynal/-ynenitrile/-ynamide/…）保留炔位次，
与对应的 prop-2-enal/prop-2-enenitrile 一致；开链烃 propyne 仍省略位次。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # alkynols (OH principal; …-n-yn-m-ol)
    ("C#CCCO", "but-3-yn-1-ol", "丁-3-炔-1-醇"),
    ("C#CCO", "prop-2-yn-1-ol", "丙-2-炔-1-醇"),
    ("CC#CCO", "but-2-yn-1-ol", "丁-2-炔-1-醇"),
    ("C#CC(O)C", "but-3-yn-2-ol", "丁-3-炔-2-醇"),
    # alkynals
    ("C#CCC=O", "but-3-ynal", "丁-3-炔醛"),
    ("C#CC=O", "prop-2-ynal", "丙-2-炔醛"),
    # alkynenitriles
    ("C#CCC#N", "but-3-ynenitrile", "丁-3-炔腈"),
    ("C#CC#N", "prop-2-ynenitrile", "丙-2-炔腈"),
    # alkynamides
    ("C#CCC(=O)N", "but-3-ynamide", "丁-3-炔酰胺"),
    ("C#CC(=O)N", "prop-2-ynamide", "丙-2-炔酰胺"),
    # negatives: ene alcohol, sat alcohol, alkynoic acid, free alkyne
    ("C/C=C/CO", "(2E)-but-2-en-1-ol", "(2E)-丁-2-烯-1-醇"),
    ("CCCCO", "butan-1-ol", "丁-1-醇"),
    ("C#CCC(=O)O", "but-3-ynoic acid", "丁-3-炔酸"),
    ("CC#C", "propyne", "丙炔"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkyne_fgs(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
