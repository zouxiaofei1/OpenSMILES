# IUPAC: P-29.3 / P-44
# Layer: L2,L3
"""Unsubstituted pyridin-n-yl side chains on open FG parents.

Chain alcohol/amine/nitrile/ketone claim pyridin-2/3/4-yl; pyridine as
parent with alkyl sides unchanged.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    (
        "OCCc1ccncc1",
        "2-(pyridin-4-yl)ethanol",
        "2-(吡啶-4-基)乙醇",
    ),
    (
        "NCCc1ccncc1",
        "2-(pyridin-4-yl)ethanamine",
        "2-(吡啶-4-基)乙胺",
    ),
    (
        "N#CCc1ccncc1",
        "2-(pyridin-4-yl)acetonitrile",
        "2-(吡啶-4-基)乙腈",
    ),
    (
        "N#CCc1ccccn1",
        "2-(pyridin-2-yl)acetonitrile",
        "2-(吡啶-2-基)乙腈",
    ),
    (
        "O=C(C)c1ccncc1",
        "1-(pyridin-4-yl)ethanone",
        "1-(吡啶-4-基)乙酮",
    ),
    (
        "OCCc1ccccn1",
        "2-(pyridin-2-yl)ethanol",
        "2-(吡啶-2-基)乙醇",
    ),
    # negatives: pyridine parent, phenyl ethanol
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Cc1ccncc1", "4-methylpyridine", "4-甲基吡啶"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine", "4-苯基吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinyl_side(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
