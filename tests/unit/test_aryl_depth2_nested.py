# IUPAC: P-29.3
# Layer: L2,L3
"""Depth-2: propoxy/butoxy leaves; nested unsub Ph on Ph arm; phenol vs chain amine.

Parent stays chain alcohol/amine; nested Ph is a recursive leaf (not Ph-on-Ph arm rejection).
Bare arene propoxy and biphenyl-as-benzene must stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # propoxy / butoxy leaves on phenylethanol
    ("CCCOc1ccc(CCO)cc1", "2-(4-propoxyphenyl)ethanol", "2-(4-丙氧基苯基)乙醇"),
    ("CCCCOc1ccc(CCO)cc1", "2-(4-butoxyphenyl)ethanol", "2-(4-丁氧基苯基)乙醇"),
    # nested unsubstituted phenyl on arm Ph
    (
        "c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-(4-phenylphenyl)ethanol",
        "2-(4-苯基苯基)乙醇",
    ),
    (
        "c1ccc(-c2cccc(CCO)c2)cc1",
        "2-(3-phenylphenyl)ethanol",
        "2-(3-苯基苯基)乙醇",
    ),
    # phenolic OH must not steal aliphatic amine parent
    (
        "Oc1ccc(CCN)cc1",
        "2-(4-hydroxyphenyl)ethanamine",
        "2-(4-羟基苯基)乙胺",
    ),
    # negatives: arene propoxy parent; biphenyl as benzene+phenyl; prior depth-2
    ("CCCOc1ccccc1", "propoxybenzene", "丙氧基苯"),
    ("CCCOc1ccc(CC)cc1", "1-ethyl-4-propoxybenzene", "1-乙基-4-丙氧基苯"),
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("Oc1ccc(CCO)cc1", "2-(4-hydroxyphenyl)ethanol", "2-(4-羟基苯基)乙醇"),
]
