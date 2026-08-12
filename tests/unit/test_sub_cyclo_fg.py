# IUPAC: P-63.1.1 / P-64.2.1 / P-62.2.1
# Layer: L2,L4,L5
"""Saturated monocyclic mono-alcohol/ketone/amine with simple ring subs.

Parent = unfused C3–C10 sat carbocycle with exactly one ring-carbon FG
(OH / ketone carbonyl / primary amine). Ring may carry halo and/or claimable
n-alkyl (C1–C4) sides. With any ring sub, keep FG locant 1
(e.g. 2-methylcyclohexan-1-ol); unsubstituted still omits locant
(cyclohexanol / cyclohexanone / cyclohexanamine).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: monoalkyl cycloalcohols
    ("CC1CCCCC1O", "2-methylcyclohexan-1-ol", "2-甲基环己-1-醇"),
    ("CCC1CCCCC1O", "2-ethylcyclohexan-1-ol", "2-乙基环己-1-醇"),
    ("CC1CCC(O)CC1", "4-methylcyclohexan-1-ol", "4-甲基环己-1-醇"),
    ("CC1CCCC1O", "2-methylcyclopentan-1-ol", "2-甲基环戊-1-醇"),
    # positive: monohalo cycloalcohol
    ("ClC1CCCCC1O", "2-chlorocyclohexan-1-ol", "2-氯环己-1-醇"),
    ("BrC1CCCCC1O", "2-bromocyclohexan-1-ol", "2-溴环己-1-醇"),
    # positive: monoalkyl / monohalo cycloamine
    ("CC1CCCCC1N", "2-methylcyclohexan-1-amine", "2-甲基环己-1-胺"),
    ("CC1CCC(N)CC1", "4-methylcyclohexan-1-amine", "4-甲基环己-1-胺"),
    ("ClC1CCCCC1N", "2-chlorocyclohexan-1-amine", "2-氯环己-1-胺"),
    # negative: unsubstituted keep omit-locant forms; open-chain / cycloalkane intact
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("NC1CCCCC1", "cyclohexanamine", "环己胺"),
    ("CCO", "ethanol", "乙醇"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sub_cyclo_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methyl_cyclohexanol_not_ethane() -> None:
    """Ring monoalcohol + methyl must not collapse to ethane."""
    r = SMILESNNamer().name("CC1CCCCC1O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylcyclohexan-1-ol"
    assert en != "ethane"
    assert "hexanol" not in en or en.startswith("2-methylcyclo")
