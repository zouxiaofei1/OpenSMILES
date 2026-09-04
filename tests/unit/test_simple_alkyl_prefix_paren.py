# IUPAC: P-16.5.1.1
# Layer: L3
"""Linear n-alkyl as substituent prefix: NO parentheses when simple & locant-free.

Regression for C5+ n-alkyl (pentyl/hexyl/heptyl) that miss the anchored table
(which only registers up to n-butyl) and fall to the recursive cut→radical-yl
backend. That backend wrongly marked every non-phenyl/non-alkoxy radical-yl
paren=True, producing "(pentyl)benzene". Simple unsubstituted alkyl stems get
no parentheses (P-16.5.1.1: only composite/complex prefixes are enclosed).

Negatives: branched (2,2-dimethylpropyl / 2-methylpropyl) and halo-alkyl
(3-bromopropyl) stems keep their parentheses — they are composite prefixes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: locant-free simple n-alkyl — no parentheses
    ("CCCCCc1ccccc1", "pentylbenzene", "戊基苯"),
    ("CCCCCCc1ccccc1", "hexylbenzene", "己基苯"),
    ("CCCCCCCc1ccccc1", "heptylbenzene", "庚基苯"),
    ("C1(CCCCC)CCCCC1", "pentylcyclohexane", "戊基环己烷"),
    # negative: composite (locant-carrying / halo) prefixes keep parentheses
    ("CC(C)(C)Cc1ccccc1", "(2,2-dimethylpropyl)benzene", "(2,2-二甲基丙基)苯"),
    ("CC(C)Cc1ccccc1", "(2-methylpropyl)benzene", "(2-甲基丙基)苯"),
    ("BrCCCc1ccccc1", "(3-bromopropyl)benzene", "(3-溴丙基)苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_alkyl_prefix_paren(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
