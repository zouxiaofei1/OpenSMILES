# IUPAC: P-29.3 / P-29.6
# Layer: L2,L3
"""Branched alkyl side chains: isobutyl, sec-butyl, neopentyl, isopentyl.

Retained prefixes (P-29.3 / P-29.6) on ring parents via shared L2 `_side_atoms`
and L3 naming. Linear C1–C4, isopropyl, and tert-butyl must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer
from namepredict.layer3.substituent_extractor import alkyl_alpha_key

CASES = [
    # positive: mono branched alkylbenzenes
    ("CC(C)Cc1ccccc1", "isobutylbenzene", "异丁基苯"),
    ("CCC(C)c1ccccc1", "sec-butylbenzene", "仲丁基苯"),
    ("CC(C)(C)Cc1ccccc1", "neopentylbenzene", "新戊基苯"),
    ("CC(C)CCc1ccccc1", "isopentylbenzene", "异戊基苯"),
    # multi-sub arene + phenol + cycloalkane (shared _side_atoms)
    ("Cc1ccc(CC(C)C)cc1", "1-isobutyl-4-methylbenzene", "1-异丁基-4-甲基苯"),
    ("Oc1ccc(CC(C)C)cc1", "4-isobutylphenol", "4-异丁基苯酚"),
    ("CC(C)CC1CCCCC1", "isobutylcyclohexane", "异丁基环己烷"),
    # negative: keep existing correct names
    ("CC(C)c1ccccc1", "isopropylbenzene", "异丙基苯"),
    ("CC(C)(C)c1ccccc1", "tert-butylbenzene", "叔丁基苯"),
    ("CCCCc1ccccc1", "butylbenzene", "丁基苯"),
    ("CC(C)CC", "2-methylbutane", "2-甲基丁烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_branched_alkyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("CC(C)Cc1ccccc1", "2-methylnonane"),
        ("CCC(C)c1ccccc1", "3-methylnonane"),
        ("CC(C)(C)Cc1ccccc1", "2,2-dimethylnonane"),
        ("CC(C)CCc1ccccc1", "2-methyldecane"),
    ],
)
def test_not_chain_alkane_name(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) != normalize_en(bad)


@pytest.mark.parametrize(
    "stem,key",
    [
        ("tert-butyl", "butyl"),
        ("sec-butyl", "butyl"),
        ("isobutyl", "isobutyl"),
        ("isopropyl", "isopropyl"),
        ("methyl", "methyl"),
    ],
)
def test_alkyl_alpha_key_strips_italic_prefixes(stem: str, key: str) -> None:
    assert alkyl_alpha_key(stem) == key
