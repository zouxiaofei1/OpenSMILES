# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained 1H-indazole scope (P-22.2.1) — negative guard only.

Positive indazole cases are not covered in this file. The retained case asserts
indole must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: must not regress indole / benzofuran / pyrazole / quinoline
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
