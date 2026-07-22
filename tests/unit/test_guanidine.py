# IUPAC: P-66.4.1.2.1
# Layer: L1,L2,L5
"""Guanidine functional parent H2N–C(=NH)–NH2 (retained name guanidine).

Scope (first cut):
- unsubstituted guanidine
- mono N-aryl guanidine (unfused Ph, Me/halo ≤2) → 1-phenylguanidine style
- N-arylsulfonyl guanidine Ar–SO2–NH–C(=NH)NH2 → 1-(…sulfonyl)guanidine dual
Must not regress methanamine, aniline, urea, benzenesulfonamide, sulfonyl chloride.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("N=C(N)N", "guanidine", "胍"),
    # positive: mono N-aryl
    ("N=C(N)Nc1ccccc1", "1-phenylguanidine", "1-苯基胍"),
    # positive: N-arylsulfonyl (benchmark dual)
    (
        "ClC1=C(C=CC=C1)S(=O)(=O)NC(=N)N",
        "1-(2-chlorophenylsulfonyl)guanidine",
        "1-(2-氯苯基磺酰基)胍",
    ),
    # negative: primary amine / aniline / urea / sulfonamide / sulfonyl chloride
    ("CN", "methanamine", "甲胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
    ("NC(=O)N", "urea", "脲"),
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride", "甲磺酰氯"),
]
