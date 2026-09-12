# IUPAC: P-29.3.2
# Layer: L2,L3,L5
"""2-methylbutan-2-yl (tert-pentyl / 1,1-dimethylpropyl) branched alkyl.

IUPAC P-29.3.2 systematic branched alkyl prefix. Topology:
  parent–C(CH3)(CH3)–CH2–CH3
Shared via L2 `_side_atoms` / L3 `_BRANCH_CHECKS`. L5 wraps leading-locant
stems in parentheses. tert-butyl and other retained branches must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.tools.re import alkyl_alpha_key
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: mono 2-methylbutan-2-yl cycloalkane / arene
    (
        "C1CC(C(C)(CC)C)CCC1",
        "(2-methylbutan-2-yl)cyclohexane",
        "(2-甲基丁-2-基)环己烷",
    ),
    (
        "CCC(C)(C)C1CCCCC1",
        "(2-methylbutan-2-yl)cyclohexane",
        "(2-甲基丁-2-基)环己烷",
    ),
    (
        "CC(C)(CC)c1ccccc1",
        "(2-methylbutan-2-yl)benzene",
        "(2-甲基丁-2-基)苯",
    ),
    (
        "CCC(C)(C)c1ccccc1",
        "(2-methylbutan-2-yl)benzene",
        "(2-甲基丁-2-基)苯",
    ),
    # negative: near-miss retained / linear must not break
    ("CC(C)(C)C1CCCCC1", "tert-butylcyclohexane", "叔丁基环己烷"),
    ("CC(C)(C)c1ccccc1", "tert-butylbenzene", "叔丁基苯"),
    ("CCC", "propane", "丙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_tert_pentyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("C1CC(C(C)(CC)C)CCC1", "cyclohexane"),
        ("CCC(C)(C)C1CCCCC1", "cyclohexane"),
        ("CC(C)(CC)c1ccccc1", "3,3-dimethylnonane"),
    ],
)
def test_not_silent_or_chain(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) != normalize_en(bad)


@pytest.mark.parametrize(
    "stem,key",
    [
        ("2-methylbutan-2-yl", "methylbutan-2-yl"),
        ("tert-butyl", "butyl"),
        ("sec-butyl", "butyl"),
    ],
)
def test_alkyl_alpha_key_tert_pentyl(stem: str, key: str) -> None:
    assert alkyl_alpha_key(stem) == key
