# IUPAC: P-22.2.1
# Layer: L2
"""Simple pyridine retains alkoxy/nitro ring substituents like benzene (P-22.2.1 / P-61.5 / P-63.2.2)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: halo + alkoxy
    ("ClC1=NC=CC=C1OC", "2-chloro-3-methoxypyridine", "2-氯-3-甲氧基吡啶"),
    # positive: multi-halo + ethoxy
    ("BrC1=NC(=CC=C1OCC)Br", "2,6-dibromo-3-ethoxypyridine", "2,6-二溴-3-乙氧基吡啶"),
    # positive: alkoxy + alkyl + nitro
    ("COC1=NC=C(C=C1C)[N+](=O)[O-]", "2-methoxy-3-methyl-5-nitropyridine", "2-甲氧基-3-甲基-5-硝基吡啶"),
    # positive probes: mono-alkoxy / mono-nitro
    ("COc1ccccn1", "2-methoxypyridine", "2-甲氧基吡啶"),
    ("[O-][N+](=O)c1ccncc1", "4-nitropyridine", "4-硝基吡啶"),
    # negative near-miss: existing simple pyridine / benzene stay correct
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Clc1ccncc1", "4-chloropyridine", "4-氯吡啶"),
    ("Cc1ccncc1", "4-methylpyridine", "4-甲基吡啶"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridine_alkoxy_nitro(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
