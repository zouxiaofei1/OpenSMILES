# IUPAC: P-22.2.1
# Layer: L2
"""Simple pyridine retains alkoxy/nitro ring substituents like benzene (P-22.2.1 / P-61.5 / P-63.2.2)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative near-miss: existing simple pyridine / benzene stay correct
    ("c1ccncc1", "pyridine", "吡啶"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridine_alkoxy_nitro(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
