# IUPAC: P-66.1.1 / P-22.2.1 / P-65.1.1.1
# Layer: L2,L3,L4,L5
"""Retained pyridinecarboxamide parent (pyridine-n-carboxamide / 吡啶-n-甲酰胺).

Unfused pyridine + single exocyclic –CONH2 (P-66.1.1 / P-22.2.1). N=1;
carboxamide attach locant; ring simple prefixes aligned with pyridinecarboxylic;
N is H or simple N-C1–C4 / N,N-dialkyl. Prefer system name over isonicotinamide.
Block (pyridin-n-yl)formamide collapse.

已知局限：pyridine-carbonitrile 的碳腈后缀（pyridine-4-carbonitrile）尚未实现。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted pyridine-2/3/4-carboxamide
    ("NC(=O)c1ccccn1", "pyridine-2-carboxamide", "吡啶-2-甲酰胺"),
    ("NC(=O)c1cccnc1", "pyridine-3-carboxamide", "吡啶-3-甲酰胺"),
    ("NC(=O)c1ccncc1", "pyridine-4-carboxamide", "吡啶-4-甲酰胺"),
    ("c1ccncc1C(=O)N", "pyridine-3-carboxamide", "吡啶-3-甲酰胺"),
    # positive: simple N-alkyl
    ("c1ccncc1C(=O)NC", "N-methylpyridine-3-carboxamide", "N-甲基吡啶-3-甲酰胺"),
    ("C(C)NC(=O)C1=CC=NC=C1", "N-ethylpyridine-4-carboxamide", "N-乙基吡啶-4-甲酰胺"),
    ("c1ccncc1C(=O)N(C)C", "N,N-dimethylpyridine-3-carboxamide", "N,N-二甲基吡啶-3-甲酰胺"),
    # positive: ring prefix (N=1; CONH2 + halo lowest set, P-14.4(c) 主官能团最低位次)
    ("O=C(N)c1ccc(Cl)nc1", "6-chloropyridine-3-carboxamide", "6-氯吡啶-3-甲酰胺"),
    ("O=C(N)c1cc(C)ncc1", "2-methylpyridine-4-carboxamide", "2-甲基吡啶-4-甲酰胺"),
    # negative: near-miss — must not steal benzamide / acetamide / acid / amine
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("O=C(O)c1ccncc1", "pyridine-4-carboxylic acid", "吡啶-4-甲酸"),
    ("Nc1ccncc1", "pyridin-4-amine", "吡啶-4-胺"),
    ("O=C(N)C1CCCCC1", "cyclohexanecarboxamide", "环己烷甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinecarboxamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
