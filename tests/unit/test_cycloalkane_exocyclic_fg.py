# IUPAC: P-66.6.1 / P-66.5.1 / P-66.1.1 / P-65.6 / P-65.5
# Layer: L2,L4,L5
"""Monocyclic cycloalkane with one exocyclic carbonyl-class FG (P-65/P-66).

Parent = unfused saturated monocarbocycle C3–C10 with exactly one ring-attached
exocyclic FG of type: carbaldehyde / carbonitrile / carboxamide /
carboxylate (alkyl ester) / carbonyl halide.
Ring may carry 0–1 simple sub: monohalo or mono C1–C4 n-alkyl; FG attach = 1.
IUPAC PIN stems: cyclohexanecarbaldehyde / carbonitrile / carboxamide /
carboxylate / carbonyl chloride (aligned with cycloalkanecarboxylic acid).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: open-chain formic-style must stay; arene / cyclo acid unchanged
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC#N", "acetonitrile", "乙腈"),
    ("CC(N)=O", "acetamide", "乙酰胺"),
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkane_exocyclic_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
