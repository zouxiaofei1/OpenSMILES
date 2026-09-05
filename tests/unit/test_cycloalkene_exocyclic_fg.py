# IUPAC: P-31.1.2 / P-65.2.2.1 / P-66.6.1.1.3
# Layer: L4,L5
"""Monocyclic cycloalkene with one exocyclic carbonyl-class FG (acid/ester/aldehyde).

Parent = unfused unsaturated monocarbocycle carrying one exocyclic COOH / COOR / CHO
(ring C=C must survive in the stem: single ene locant-1 is elided in EN but kept in ZH
cyclohexene vs 环己-1-烯; non-1 / polyene locants explicit). FG locant becomes explicit
once the ring is unsaturated because numbering is no longer unique.
Saturated cycloalkanecarboxylic acid / open-chain / arene must not change.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: single endocyclic C=C, ene locant 1 (EN elided), FG@1 explicit
    ("O=C(O)C1=CC(=O)[C@@H](O)[C@H](O)C1",
     "(4S,5R)-4,5-dihydroxy-3-oxocyclohexene-1-carboxylic acid",
     "(4S,5R)-4,5-二羟基-3-氧代环己-1-烯-1-羧酸"),
    # carboxylate anion (去质子酸根)
    ("O=C([O-])C1=CC(=O)[C@@H](O)[C@H](O)C1",
     "(4S,5R)-4,5-dihydroxy-3-oxocyclohexene-1-carboxylate",
     "(4S,5R)-4,5-二羟基-3-氧代环己-1-烯-1-羧酸根"),
    # positive: diene, polyene locants explicit
    ("N[C@@H]1C(C(=O)O)=CC=C[C@@H]1O",
     "(5S,6R)-6-amino-5-hydroxycyclohexa-1,3-diene-1-carboxylic acid",
     "(5S,6R)-6-氨基-5-羟基环己-1,3-二烯-1-羧酸"),
    ("N[C@H]1C(C(=O)O)=CC=C[C@@H]1O",
     "(5S,6S)-6-amino-5-hydroxycyclohexa-1,3-diene-1-carboxylic acid",
     "(5S,6S)-6-氨基-5-羟基环己-1,3-二烯-1-羧酸"),
    ("CC1=CC=CC(O)C1(O)C(=O)O",
     "1,6-dihydroxy-2-methylcyclohexa-2,4-diene-1-carboxylic acid",
     "1,6-二羟基-2-甲基环己-2,4-二烯-1-羧酸"),
    # positive: ester / aldehyde, ene locant non-1 / 1
    ("CCOC(=O)[C@@]1(c2ccccc2)CCC=C[C@@H]1N(C)C",
     "ethyl (1R,2S)-2-(dimethylamino)-1-phenylcyclohex-3-ene-1-carboxylate",
     None),  # zh 二甲氨基 vs 二甲基氨基 为独立胺名规则，不在 R1 范围
    ("C=C(C)C1CC=C(C=O)CC1",
     "4-prop-1-en-2-ylcyclohexene-1-carbaldehyde",
     "4-丙-1-烯-2-基环己-1-烯-1-甲醛"),
    # negative: saturated / open-chain / arene stay intact
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷羧酸"),
    ("OC(=O)C1CCCC1", "cyclopentanecarboxylic acid", "环戊烷羧酸"),
    ("OC(=O)CCCCCC", "heptanoic acid", "庚酸"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkene_exocyclic_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
