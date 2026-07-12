# IUPAC: P-66.5.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated mononitriles (alkenenitriles).

Nitrile is the principal characteristic group (C≡N carbon = locant 1); one
non-aromatic C=C is expressed as -n-enenitrile / -n-烯腈 with the lower
double-bond carbon locant. No (E)/(Z) stereodescriptors this cycle.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenenitriles (no E/Z)
    ("C=CC#N", "prop-2-enenitrile", "丙-2-烯腈"),
    ("CC=CC#N", "but-2-enenitrile", "丁-2-烯腈"),
    ("C=CCC#N", "but-3-enenitrile", "丁-3-烯腈"),
    ("CCC=CC#N", "pent-2-enenitrile", "戊-2-烯腈"),
    ("C=CCCC#N", "pent-4-enenitrile", "戊-4-烯腈"),
    # negative: saturated nitriles, alkenoic acid, alkenal must not become alkenenitriles
    ("CC#N", "acetonitrile", "乙腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenenitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
