# IUPAC: P-14.2.1
# Layer: L5
"""Multiplicative prefixes penta–deca for identical substituents (P-14.2.1 / P-16.3).

Affiliated: P-61.3.1 halogen stems; existing di/tri/tetra must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: count ≥5 needs penta/hexa/… multiplicative prefixes
    # zh=None when gold has no Chinese (merged_benchmark empty chinese_name)
    ("ClC(Cl)C(Cl)(Cl)Cl", "1,1,1,2,2-pentachloroethane", None),
    ("ClC(Cl)(Cl)C(Cl)Cl", "1,1,1,2,2-pentachloroethane", None),
    ("ClC(Cl)(Cl)C(Cl)(Cl)Cl", "1,1,1,2,2,2-hexachloroethane", None),
    ("FC(F)(F)C(F)(F)F", "1,1,1,2,2,2-hexafluoroethane", None),
    ("FC(F)(F)C(F)(F)C(F)F", "1,1,1,2,2,3,3-heptafluoropropane", None),
    ("FC(F)(F)C(F)(F)C(F)(F)F", "1,1,1,2,2,3,3,3-octafluoropropane", None),
    # negative: di/tetra already correct — must not regress
    ("ClC(Cl)(Cl)Cl", "tetrachloromethane", "四氯甲烷"),
    ("ClCCCl", "1,2-dichloroethane", "1,2-二氯乙烷"),
    ("BrCBr", "dibromomethane", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mult_prefix(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
