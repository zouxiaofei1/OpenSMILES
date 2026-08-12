# IUPAC: P-66.1.1.1.3
# Layer: L2,L3,L5
"""Open-chain N-phenyl amide scope (P-66.1.1.1.3) — negative guard only.

Positive N-phenyl amide cases are not covered in this file. The retained cases
assert primary amide and C-phenyl benzamide are not named as N-phenyl amides.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: N-alkyl / primary stay correct; C-phenyl benzamide not N-phenyl
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccc(C(=O)N)cc1", "benzamide", "苯甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_phenyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
