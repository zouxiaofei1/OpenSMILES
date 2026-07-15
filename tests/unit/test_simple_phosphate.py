# IUPAC: P-67
# Layer: L1,L2,L5
"""Simple monoalkyl dihydrogen phosphate and alkylphosphonic acid (P-67).

Minimal first cut: open-chain monoalkyl phosphate monoester and mono C–P
phosphonic acid; no nucleotides / polyphosphates / salt class names.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: monoalkyl dihydrogen phosphates
    ("COP(=O)(O)O", "methyl dihydrogen phosphate", "磷酸甲酯"),
    ("CCOP(=O)(O)O", "ethyl dihydrogen phosphate", "磷酸乙酯"),
    ("CCCCCCCCCCCCCOP(=O)(O)O", "tridecyl dihydrogen phosphate", "磷酸十三烷基酯"),
    # positive: simple alkylphosphonic acid
    ("CP(=O)(O)O", "methylphosphonic acid", "甲基膦酸"),
    # negative: must not misclassify carboxylic acid / alcohol
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_phosphate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
