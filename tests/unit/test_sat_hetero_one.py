# IUPAC: P-65.6.3.5.1 / P-66.1.5.1
# Layer: L2,L4,L5
"""Sat-monohetero lactone / lactam scope (P-65.6.3.5.1) — negative guard only.

Positive lactone/lactam cases are not covered in this file. The retained cases
assert open-chain ester / amide keep their open-chain names.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: open-chain ester / amide keep open-chain names
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_one(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
