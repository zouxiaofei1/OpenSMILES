# IUPAC: P-66.1.1.1.3
# Layer: L1,L2,L3,L5
"""N-羟基(异羟肟酸型)酰胺 scope — 正例 + 近邻负例。

N-羟基酰胺/异羟肟酸: 羰基 C(=O) 上的酰胺 N 除羰基碳外还带一个不再连碳的
羟基氧(N-OH),如 acetohydroxamic acid 按 N-hydroxy-…amide 表达。此前 L1
`_amide_n_info` 只允许酰胺 N 上的取代基为 C/H,N-OH 使酰胺漏检、羰基碳被误判
成醛 → 输出 1-[…]acetaldehyde(P-66.1.1.1.3 N-前缀命名)。

负例: N-烷基酰胺、伯酰胺、酮、酸、醛不得被误伤。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: N-hydroxy amide family (gold chebi-1411)
    ("CC(=O)N(O)CCCN", "N-(3-aminopropyl)-N-hydroxyacetamide", "N-(3-氨基丙基)-N-羟基乙酰胺"),
    ("CC(=O)NO", "N-hydroxyacetamide", "N-羟基乙酰胺"),
    ("CC(=O)N(O)C", "N-hydroxy-N-methylacetamide", "N-羟基-N-甲基乙酰胺"),
    ("CC(=O)N(C)O", "N-hydroxy-N-methylacetamide", "N-羟基-N-甲基乙酰胺"),
    ("CC(=O)N(O)O", "N,N-dihydroxyacetamide", "N,N-二羟基乙酰胺"),
    # negative: N-alkyl amide / primary amide / ketone / acid / aldehyde stay correct
    ("CC(=O)NCCCN", "N-(3-aminopropyl)acetamide", "N-(3-氨基丙基)乙酰胺"),
    ("CC(=O)NC", "N-methylacetamide", "N-甲基乙酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CC(=O)CCC", "pentan-2-one", "戊-2-酮"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_hydroxy_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
