# IUPAC: P-63.2.2
# Layer: L2, L5
"""Branched dialkyl ethers: isopropyl / hexafluoroisopropyl functional class.

Linear C1–C4 arms keep alkoxyalkane / retained sym ether. Branched specials
use functional-class names (… methyl ether / …基甲醚).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # hexafluoroisopropyl methyl ether (benchmark tiers-15987)
    (
        "COC(C(F)(F)F)C(F)(F)F",
        "hexafluoroisopropyl methyl ether",
        "六氟异丙基甲醚",
    ),
    (
        "FC(F)(F)C(OC)C(F)(F)F",
        "hexafluoroisopropyl methyl ether",
        "六氟异丙基甲醚",
    ),
    # unsubstituted isopropyl methyl ether
    ("COC(C)C", "isopropyl methyl ether", "异丙基甲醚"),
    # linear regressions
    ("COC", "dimethyl ether", "二甲基醚"),
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCOCC", "diethyl ether", "二乙基醚"),
    # negatives
    ("CCO", "ethanol", "乙醇"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
]
