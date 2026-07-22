# IUPAC: P-29.3 / P-25
# Layer: L2,L3
"""Unsubstituted naphthalen-n-yl side chains on open FG parents."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    (
        "OCCc1ccc2ccccc2c1",
        "2-(naphthalen-2-yl)ethanol",
        "2-(萘-2-基)乙醇",
    ),
    (
        "OCCc1cccc2ccccc12",
        "2-(naphthalen-1-yl)ethanol",
        "2-(萘-1-基)乙醇",
    ),
    (
        "N#CCc1ccc2ccccc2c1",
        "2-(naphthalen-2-yl)acetonitrile",
        "2-(萘-2-基)乙腈",
    ),
    (
        "NCCc1cccc2ccccc12",
        "2-(naphthalen-1-yl)ethanamine",
        "2-(萘-1-基)乙胺",
    ),
    (
        "O=C(C)c1ccc2ccccc2c1",
        "1-(naphthalen-2-yl)ethanone",
        "1-(萘-2-基)乙酮",
    ),
    # negatives
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("Cc1ccc2ccccc2c1", "2-methylnaphthalene", "2-甲基萘"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("OCCc1ccncc1", "2-(pyridin-4-yl)ethanol", "2-(吡啶-4-基)乙醇"),
]
