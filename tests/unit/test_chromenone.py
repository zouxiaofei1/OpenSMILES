# IUPAC: P-25 / P-22.2.1 / P-65.6.3
# Layer: L2,L3,L4,L5
"""Retained parent chromen-2-one (coumarin lactone; 2H-1-benzopyran-2-one style).

Fused 6+6: benzene + α-pyrone lactone (9C + ring O + exocyclic =O at C2).
O=1, carbonyl C=2; EN stem chromen-2-one; ZH 香豆素. Simple ring prefixes
(halo / Me / alkoxy / OH / amino / nitro / prenyl / phenyl) only.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted retained parent
    ("O=c1ccc2ccccc2o1", "chromen-2-one", "香豆素"),
    # positive: simple ring prefixes
    ("COc1ccc2ccc(=O)oc2c1", "7-methoxychromen-2-one", "7-甲氧基香豆素"),
    ("Cc1ccc2ccc(=O)oc2c1", "7-methylchromen-2-one", "7-甲基香豆素"),
    ("Clc1ccc2ccc(=O)oc2c1", "7-chlorochromen-2-one", "7-氯香豆素"),
    # positive: anchor (prenyl + methoxy)
    (
        "COc1cc2oc(=O)ccc2cc1CC=C(C)C",
        "7-methoxy-6-(3-methylbut-2-enyl)chromen-2-one",
        "7-甲氧基-6-(3-甲基丁-2-烯基)香豆素",
    ),
    # positive: optional phenyl (omit 2H- to match project stem)
    ("c1ccc(cc1)c1cc(=O)oc2ccccc12", "4-phenylchromen-2-one", "4-苯基香豆素"),
    # negative: must not become chromen-2-one
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
]
