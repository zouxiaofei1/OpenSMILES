# IUPAC: P-66.6.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monoaldehydes (alkenals).

Aldehyde is the principal characteristic group (CHO = locant 1); one non-aromatic
C=C is expressed as -n-enal / -n-烯醛 with the lower double-bond carbon locant.
No (E)/(Z) stereodescriptors this cycle.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenals (no E/Z)
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("C=CC=O", "prop-2-enal", "丙-2-烯醛"),
    ("CCC=CC=O", "pent-2-enal", "戊-2-烯醛"),
    ("CCCC=CC=O", "hex-2-enal", "己-2-烯醛"),
    ("C=CCC=O", "but-3-enal", "丁-3-烯醛"),
    # negative: saturated aldehydes and alkenoic acid must not become alkenals
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCCC=O", "butanal", "丁醛"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenal(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
