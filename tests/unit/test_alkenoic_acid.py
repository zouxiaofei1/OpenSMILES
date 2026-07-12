# IUPAC: P-65.1.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monocarboxylic acids (alkenoic acids).

Carboxyl is principal characteristic group (COOH = locant 1); one non-aromatic
C=C is expressed as -n-enoic / -n-烯酸 with the lower double-bond carbon locant.
No (E)/(Z) stereodescriptors this cycle.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenoic acids (no E/Z)
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    ("C=CCC(=O)O", "but-3-enoic acid", "丁-3-烯酸"),
    ("CC=CC(=O)O", "but-2-enoic acid", "丁-2-烯酸"),
    ("CC(C)=CC(=O)O", "3-methylbut-2-enoic acid", "3-甲基丁-2-烯酸"),
    ("CC/C=C/CC(=O)O", "hex-3-enoic acid", "己-3-烯酸"),
    ("CCC=CC(=O)O", "pent-2-enoic acid", "戊-2-烯酸"),
    # negative: saturated acid, alkene, oxoacid must not become alkenoic acids
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("C=C", "ethene", "乙烯"),
    ("OC(=O)C(=O)C", "2-oxopropanoic acid", "2-氧代丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
