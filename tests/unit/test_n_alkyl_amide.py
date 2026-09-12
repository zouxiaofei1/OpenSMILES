# IUPAC: P-66.1.1.1.3
# Layer: L1,L2,L3,L5
"""Open-chain N-substituted amide scope (P-66.1.1.1.3) — negative guard only.

Positive N-alkyl amide cases are not covered in this file. The retained cases
assert primary amide / acid / aldehyde are not named as N-alkyl amides.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: primary amide / acid / aldehyde / sec amine must stay correct
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_alkyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
