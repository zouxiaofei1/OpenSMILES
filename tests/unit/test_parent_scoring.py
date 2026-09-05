# IUPAC: P-44.1 / P-29.3.1
# Layer: L2,L3
"""Composable parent scoring + recursive alkyl substituents.

Layer2 select_parent moves from an if-else waterfall to candidate
collection + composable scoring: a benzene (or cycloalkane) core with
fully-recognizable alkyl sides must win over the longest-chain alkane,
while principal-FG parents (e.g. benzoic acid) must still outrank the
bare ring. Layer3 recognizes tert-butyl recursively (P-29.3.1).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: ring core must beat longest-chain alkane via scoring
    ("c1ccc(cc1)C(C)(C)C", "tert-butylbenzene", "叔丁基苯"),
    ("CCc1ccccc1C", "1-ethyl-2-methylbenzene", "1-乙基-2-甲基苯"),
    ("CC(C)(C)C1CCCCC1", "tert-butylcyclohexane", "叔丁基环己烷"),
    # negative / regression: existing correct behavior must not change
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),  # acid outranks bare benzene
    ("Cc1ccc(Cl)cc1", "1-chloro-4-methylbenzene", "1-氯-4-甲基苯"),
    ("CCC(C)C", "2-methylbutane", "2-甲基丁烷"),  # plain alkane: no ring bias
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_parent_scoring(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("c1ccc(cc1)C(C)(C)C", "2,2-dimethyloctane"),
        ("CCc1ccccc1C", "nonane"),
        # unhandled complex side: must NOT collapse to a bare "benzene"
        ("CC(C)Cc1ccccc1", "benzene"),
        ("CCC(C)c1ccccc1", "benzene"),
    ],
)
def test_no_wrong_parent(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert not (r.success and normalize_en(r.en) == normalize_en(bad))
