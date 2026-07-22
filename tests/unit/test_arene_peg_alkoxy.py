# IUPAC: P-63.2.2 / P-29.3
# Layer: L2,L3,L5
"""Arene linear n-alkoxy C1–C4 and PEG tails -(OCH2CH2)k-OR on benzene.

Extends ring alkoxy beyond methoxy/ethoxy/2-methoxyethoxy so polyether
side chains stay on the arene parent (not alkane fallback → hexane).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # simple n-alkoxy C3–C4
    ("CCCOc1ccccc1", "propoxybenzene", "丙氧基苯"),
    ("CCCCOc1ccccc1", "butoxybenzene", "丁氧基苯"),
    # PEG k=1 ethoxy
    (
        "CCOCCOc1ccccc1",
        "(2-ethoxyethoxy)benzene",
        "(2-乙氧基乙氧基)苯",
    ),
    # PEG k=2: Ph-O-(CH2CH2O)2-R  (note extra C vs mistaken COCCOCOc)
    (
        "COCCOCCOc1ccccc1",
        "(2-(2-methoxyethoxy)ethoxy)benzene",
        "(2-(2-甲氧基乙氧基)乙氧基)苯",
    ),
    (
        "CCOCCOCCOc1ccccc1",
        "(2-(2-ethoxyethoxy)ethoxy)benzene",
        "(2-(2-乙氧基乙氧基)乙氧基)苯",
    ),
    # regressions
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("CCOc1ccccc1", "ethoxybenzene", None),
    (
        "COCCOc1ccccc1",
        "(2-methoxyethoxy)benzene",
        None,
    ),
    (
        "BrC1=CC(=C(C=C1)OCCOC)F",
        "4-bromo-2-fluoro-1-(2-methoxyethoxy)benzene",
        "4-溴-2-氟-1-(2-甲氧基乙氧基)苯",
    ),
    # open-chain ether unchanged
    ("CCOCC", "diethyl ether", None),
]
