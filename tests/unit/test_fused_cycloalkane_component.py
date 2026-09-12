# IUPAC: P-25.3.2.2.1
# Layer: L2,L5
"""环烷烃作稠合附加组分：饱和单环烃名删尾 'ne' 得前缀（cyclopentane → cyclopenta），
表示最大数目非累积双键的形式（P-25.3.2.2.1）；数字位次按 P-25.3.8.1 省略。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# (smiles, EN 末尾骨架, ZH 子串)
FUSED_CASES = [
    # 环戊烷稠合四氢呋喃（无取代）：附加组分 cyclopenta，母体 furan
    ("C1CC2CCOC2C1", "cyclopenta[b]furan", "环戊并[b]呋喃"),
    # 同上 5 位羟基：位次随稠环编号（桥头 3a/6a，外周 4/5/6）
    ("O[C@@H]1C[C@@H]2OCC[C@@H]2C1", "cyclopenta[b]furan-5-ol", "环戊并[b]呋喃-5-醇"),
]


@pytest.mark.parametrize("smiles,en_tail,zh_sub", FUSED_CASES)
def test_cycloalkane_as_fused_component(smiles: str, en_tail: str, zh_sub: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert r.en.endswith(en_tail), r.en
    assert zh_sub in r.zh, r.zh


# 芳环侧稠合：附加组分与母体共享的稠合键为芳香键，碳环组分仍须被识别
AROMATIC_FUSED_CASES = [
    # 环戊烷并吡啶（共享稠合键位于吡啶 2,3 侧）
    ("c1cnc2c(c1)CCC2", "cyclopenta[b]pyridine", "环戊并[b]吡啶"),
    # 同上换母体杂环：吡嗪
    ("c1cnc2c(n1)CCC2", "cyclopenta[b]pyrazine", "环戊并[b]吡嗪"),
    # 环大小随附加组分变：环丁并/环庚并
    ("c1cnc2c(c1)CC2", "cyclobuta[b]pyridine", "环丁并[b]吡啶"),
    ("c1cnc2c(c1)CCCCC2", "cyclohepta[b]pyridine", "环庚并[b]吡啶"),
]


@pytest.mark.parametrize("smiles,en_tail,zh_sub", AROMATIC_FUSED_CASES)
def test_aromatic_fused_cycloalkane_component(smiles: str, en_tail: str, zh_sub: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert r.en.endswith(en_tail), r.en
    assert zh_sub in r.zh, r.zh


# 近邻负例：已注册稠环与单环烃不得被本轮组分路径改写
NEIGHBORS = [
    ("c1ccc2occc2c1", "1-benzofuran", "1-苯并呋喃"),
    ("C1CCC2CCCCC2C1", "decahydronaphthalene", "十氢萘"),
    ("C1CCc2ccccc2C1", "1,2,3,4-tetrahydronaphthalene", "1,2,3,4-四氢萘"),
    ("C1CCCC1", "cyclopentane", "环戊烷"),
    ("OC1CCCC1", "cyclopentanol", "环戊醇"),
    # 六元碳环稠合：整骨架有保留名模板（喹啉），不得改走环烷烃附加组分路径
    ("c1cnc2c(c1)CCCC2", "5,6,7,8-tetrahydroquinoline", "5,6,7,8-四氢喹啉"),
]


@pytest.mark.parametrize("smiles,en,zh", NEIGHBORS)
def test_neighbors_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh
