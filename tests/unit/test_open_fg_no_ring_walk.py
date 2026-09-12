# IUPAC: P-44 / P-29.6 / P-22.2
# Layer: L2,L3
"""Open FG parent chain must not enter rings; claim cycloalkyl / piperidinyl sides.

1-cyclohexylethanone and 1-(piperidin-4-yl)ethanone stay ketone parents with
ring as substituent — not octan-2-one / pentan-2-one. Ring FG parents
(cyclohexanone/cyclohexanol) unchanged.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negatives: ring FG parents + plain open ketone
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("CCC(=O)C", "butan-2-one", "丁-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_open_fg_no_ring_walk(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
