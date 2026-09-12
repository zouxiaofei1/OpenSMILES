# IUPAC: P-65.1.1
# Layer: L2,L4,L5
"""Monocyclic cycloalkanecarboxylic acids (unsub or mono ring sub).

Parent = saturated monocarbocycle with one exocyclic COOH on a ring carbon.
Allow 0–1 ring simple sub: monohalo (F/Cl/Br/I) or mono C1–C4 n-alkyl.
COOH attachment = locant 1; ring subs get lowest compatible locants.
IUPAC: cycloalkanecarboxylic acid / 环…烷甲酸.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: arene acid / open-chain acid / cycloalcohol must not change
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("OC(=O)CCCCCC", "heptanoic acid", "庚酸"),
    ("CCCCCCC(=O)O", "heptanoic acid", "庚酸"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkanecarboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
