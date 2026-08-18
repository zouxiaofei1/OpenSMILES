# IUPAC: P-65.1.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monocarboxylic alkynoic acids / alkynoates.

Carboxyl (or ester carbonyl) is principal FG (locant 1); one non-aromatic C≡C
is expressed as -n-ynoic / -n-炔酸 (or -ynoate / 炔酸…酯). First scope: mono
C≡C, no C=C, open-chain mono acid/ester.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkynoic acids
    ("C#CCC(=O)O", "but-3-ynoic acid", "丁-3-炔酸"),
    ("C#CC(=O)O", "propynoic acid", "丙炔酸"),
    ("C(CCCCCCC#C)(=O)O", "non-8-ynoic acid", "壬-8-炔酸"),
    ("C#CCC(=O)[O-]", "but-3-ynoate", "丁-3-炔酸根"),
    # positive: alkynoates (esters)
    ("C(C#C)(=O)OCC", "ethyl propynoate", "丙炔酸乙酯"),
    ("C#CC(=O)OC", "methyl propynoate", "丙炔酸甲酯"),
    # negative: saturated acid, alkenoic acid, free alkyne must not become ynoic
    ("CCCC(=O)O", "butanoic acid", "丁酸"),
    ("C/C=C/C(=O)O", "(2E)-but-2-enoic acid", "(2E)-丁-2-烯酸"),
    ("CC#C", "propyne", "丙炔"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkynoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
