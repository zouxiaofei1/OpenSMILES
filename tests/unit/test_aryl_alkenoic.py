# IUPAC: P-31.1 / P-65.1.1
# Layer: L2
"""Aryl-substituted open-chain alkenoic acids / alkenoates.

IUPAC P-31.1 / P-65.1.1 / P-65.6: the parent is the open-chain unsaturation +
principal FG (COOH/ester). An aromatic ring is a substituent, not a reason to
reject alkenoic/alkenoate selection. Gate is parent atoms not in ring, not
molecule-level has_ring.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: aryl-substituted open-chain alkenoic acids / alkenoates
    (
        r"OC(=O)/C=C/C1=CC=CC=C1",
        "(E)-3-phenylprop-2-enoic acid",
        "(E)-3-苯基丙-2-烯酸",
    ),
    (
        r"O=C(O)/C=C/c1ccc(C(F)(F)F)cc1",
        "(E)-3-[4-(trifluoromethyl)phenyl]prop-2-enoic acid",
        "(E)-3-(4-三氟甲基苯基)丙-2-烯酸",
    ),
    (
        r"COC(=O)/C=C/c1ccc(O)c(OC)c1",
        "methyl (E)-3-(4-hydroxy-3-methoxyphenyl)prop-2-enoate",
        "(E)-3-(4-羟基-3-甲氧基苯基)丙-2-烯酸甲酯",
    ),
    (
        r"COC(=O)/C=C/c1ccc(OC)c(OC)c1",
        "methyl (E)-3-(3,4-dimethoxyphenyl)prop-2-enoate",
        "(E)-3-(3,4-二甲氧基苯基)丙-2-烯酸甲酯",
    ),
    # positive: no-ring control remains alkenoic
    (r"O=C(O)/C=C/C", "(E)-but-2-enoic acid", "(E)-丁-2-烯酸"),
    # negative: saturated / diacid / ring-acid must not become enoic
    ("O=C(O)CCc1ccccc1", "3-phenylpropanoic acid", "3-苯基丙酸"),
    ("O=C(O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]
