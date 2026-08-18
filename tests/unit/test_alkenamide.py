# IUPAC: P-66.1.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monoamides (alkenamides).

Amide is the principal characteristic group (C(=O)N carbon = locant 1); one
non-aromatic C=C is expressed as -n-enamide / -n-烯酰胺 with the lower
double-bond carbon locant and (E)/(Z) when stereo is defined (P-31.1).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # P-66.1.1 / P-31.1: prop-2-enamide → retained acrylamide (IUPAC P-66.1.1.1.1)
    ("NC(=O)C=C", "acrylamide", "丙烯酰胺"),
    ("NC(=O)/C=C/C", "(2E)-but-2-enamide", "(2E)-丁-2-烯酰胺"),
    ("NC(=O)CC=C", "but-3-enamide", "丁-3-烯酰胺"),
    ("NC(=O)/C=C/c1ccccc1", "(2E)-3-phenylprop-2-enamide", "(2E)-3-苯基丙-2-烯酰胺"),
    (r"CCCCCCCC/C=C\CCCCCCCC(=O)N", "(9Z)-octadec-9-enamide", "(9Z)-十八-9-烯酰胺"),
    # negative: saturated amides must stay saturated
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("NC(=O)CCC", "butanamide", "丁酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
