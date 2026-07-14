# IUPAC: P-44 / P-29.6 / P-22.2
# Layer: L2,L3
"""Open FG parent chain must not enter rings; claim cycloalkyl / piperidinyl sides.

1-cyclohexylethanone and 1-(piperidin-4-yl)ethanone stay ketone parents with
ring as substituent — not octan-2-one / pentan-2-one. Ring FG parents
(cyclohexanone/cyclohexanol) unchanged.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # carbocycle side on open ketone
    ("O=C(C)C1CCCCC1", "1-cyclohexylethanone", "1-环己基乙酮"),
    ("CC(=O)C1CCCCC1", "1-cyclohexylethanone", "1-环己基乙酮"),
    # aza-sat side on open ketone / alcohol / amine
    ("O=C(C)C1CCNCC1", "1-(piperidin-4-yl)ethanone", "1-(哌啶-4-基)乙酮"),
    (
        "O=C(CBr)C1CCNCC1",
        "2-bromo-1-(piperidin-4-yl)ethanone",
        "2-溴-1-(哌啶-4-基)乙酮",
    ),
    ("CC(O)C1CCNCC1", "1-(piperidin-4-yl)ethanol", "1-(哌啶-4-基)乙醇"),
    ("NCC1CCNCC1", "(piperidin-4-yl)methanamine", "(哌啶-4-基)甲胺"),
    # pyrrolidine
    ("O=C(C)C1CCNC1", "1-(pyrrolidin-3-yl)ethanone", "1-(吡咯烷-3-基)乙酮"),
    # negatives: ring FG parents + plain open ketone
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("CCC(=O)C", "butan-2-one", "丁-2-酮"),
    ("O=C(C)c1ccccc1", "acetophenone", "苯乙酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_open_fg_no_ring_walk(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
