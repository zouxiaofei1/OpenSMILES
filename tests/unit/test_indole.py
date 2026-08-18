# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained 1H-indole scope (P-22.2.1) — unsubstituted + negative guard.

Only the unsubstituted indole positive case is retained here; substituted
cases are not covered. The other cases assert naphthalene / benzene / pyridine
must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
    # negative: must not regress naphthalene / imidazole / benzene / pyridine
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2[nH]ccc2c1")
    assert r.success
    assert normalize_en(r.en) == "1h-indole"
    assert "octane" not in normalize_en(r.en)
