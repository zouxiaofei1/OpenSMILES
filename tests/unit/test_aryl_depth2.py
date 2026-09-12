# IUPAC: P-29.3
# Layer: L2,L3
"""Depth-2 simple leaves on Ph/OPh/Bn arms (alkoxy / nitro / CF3).

Parent remains chain alcohol/acid/amine; substituted phenyl is the recursive arm.
Halo/methyl depth-1 and bare arene methoxy parents must stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: methoxy / ethoxy / nitro / CF3 on phenylethanol arm
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("CCOc1ccc(CCO)cc1", "2-(4-ethoxyphenyl)ethanol", "2-(4-乙氧基苯基)乙醇"),
    ("O=[N+]([O-])c1ccc(CCO)cc1", "2-(4-nitrophenyl)ethanol", "2-(4-硝基苯基)乙醇"),
    (
        "FC(F)(F)c1ccc(CCO)cc1",
        "2-[4-(trifluoromethyl)phenyl]ethanol",
        "2-(4-三氟甲基苯基)乙醇",
    ),
    # positive: acid / amine parents
    (
        "O=C(O)Cc1ccc(OC)cc1",
        "2-(4-methoxyphenyl)acetic acid",
        "2-(4-甲氧基苯基)乙酸",
    ),
    (
        "NCCc1ccc(C(F)(F)F)cc1",
        "2-[4-(trifluoromethyl)phenyl]ethanamine",
        "2-(4-三氟甲基苯基)乙胺",
    ),
    # positive: meta methoxy + benzyl-style methanol
    ("COc1cccc(CO)c1", "(3-methoxyphenyl)methanol", "(3-甲氧基苯基)甲醇"),
    # negative: depth-1 halo arm + bare arene alkoxy (not recursive)
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
    ("COc1ccc(CC)cc1", "1-ethyl-4-methoxybenzene", "1-乙基-4-甲氧基苯"),
]
