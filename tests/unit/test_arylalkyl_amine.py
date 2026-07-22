# IUPAC: P-62.2.2.1 / P-29.3
# Layer: L2,L3
"""Sec-amine: prefer aryl-bearing arm as parent; claim Ph on that arm.

When N has two equal-length alkyl arms and one carries Ph, parent is the
arylalkyl chain (not the plain ethyl). Primary 1-arylethanamine stays correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # user ex1-ish primary already ok; keep as negative
    (
        "ClC1=CC=CC=C1C(C)N",
        "1-(2-chlorophenyl)ethanamine",
        "1-(2-氯苯基)乙胺",
    ),
    (
        "FC(C1=CC=C(CCN)C=C1)(F)F",
        "2-[4-(trifluoromethyl)phenyl]ethanamine",
        "2-(4-三氟甲基苯基)乙胺",
    ),
    # user ex8 + family: sec amine + 1-arylalkyl
    (
        "C(C)NC(C)C1=CC=C(C=C1)OC",
        "N-ethyl-1-(4-methoxyphenyl)ethanamine",
        "N-乙基-1-(4-甲氧基苯基)乙胺",
    ),
    (
        "CCNC(C)c1ccccc1",
        "N-ethyl-1-phenylethanamine",
        "N-乙基-1-苯基乙胺",
    ),
    (
        "CCNCCc1ccccc1",
        "N-ethyl-2-phenylethanamine",
        "N-乙基-2-苯基乙胺",
    ),
    (
        "CNC(C)c1ccccc1",
        "N-methyl-1-phenylethanamine",
        "N-甲基-1-苯基乙胺",
    ),
    # nitriles from user set (should already work)
    (
        "C(C)C=1C=C(C=CC1)CC#N",
        "2-(3-ethylphenyl)acetonitrile",
        "2-(3-乙基苯基)乙腈",
    ),
    (
        "ClC1=C(C(=CC=C1)F)CC#N",
        "2-(2-chloro-6-fluorophenyl)acetonitrile",
        "2-(2-氯-6-氟苯基)乙腈",
    ),
    # ester: prefix inserts between alkyl and acyl (user ex7)
    (
        "COC(C(CCC1=CC=CC=C1)=O)=O",
        "methyl 2-oxo-4-phenylbutanoate",
        "2-氧代-4-苯基丁酸甲酯",
    ),
    # negatives: plain sec amine / primary
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("NCCc1ccccc1", "phenylethanamine", "苯基乙胺"),
]
