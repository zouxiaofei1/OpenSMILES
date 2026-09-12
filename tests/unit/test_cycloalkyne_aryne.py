# IUPAC: P-31.1 / P-31.2
# Layer: L2,L4,L5
"""Endocyclic alkyne carbocycle (aryne): cyclohexa-1,3-dien-5-yne.

环内三键使该环无法构成 mancude 芳香体系：RDKit 把苯炔 c1ccccc#1 的 sp 碳按 6π 一并标为
芳香，单环全芳香又匹配不到保留名时母体选择落空（no_assemblable_candidate）。环内三键
须按环烯炔表达，不得当保留芳名或未注册芳环丢弃。苯与纯烯环不得回归。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # 苯炔（C6H4）：三种等价 Kekulé 写法须同名
    ("C1=CC=CC#C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    ("C1=CC#CC=C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    ("C1#CC=CC=C1", "cyclohexa-1,3-dien-5-yne", "环己-1,3-二烯-5-炔"),
    # 带取代基的芳炔：环内三键走 carbocycle 路径，取代基照常表达
    ("CC1=CC=CC#C1", "1-methylcyclohexa-1,3-dien-5-yne", "1-甲基环己-1,3-二烯-5-炔"),
    # 非芳香环炔（RDKit 本就不标芳香）：不得受影响
    ("C1#CC=CC1", "cyclopent-1-en-3-yne", "环戊-1-烯-3-炔"),
    # 负例：真芳香环与纯烯环不得被当成烯炔
    ("C1=CC=CC=C1", "benzene", "苯"),
    ("C1=CCC=CC1", "cyclohexa-1,4-diene", "环己-1,4-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cycloalkyne_aryne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles} failed: {(r.meta or {}).get('reason')}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CC=CC#C1", "C1=CC#CC=C1", "C1#CC=CC=C1"])
def test_aryne_not_named_benzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert "benzene" not in normalize_en(r.en)
    assert "苯" != normalize_zh(r.zh)
