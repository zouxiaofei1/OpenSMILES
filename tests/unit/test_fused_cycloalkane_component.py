# IUPAC: P-25.3.2.2.1
# Layer: L2,L5
"""环烷烃作稠合附加组分：饱和单环烃名删尾 'ne' 得前缀（cyclopentane → cyclopenta），
表示最大数目非累积双键的形式（P-25.3.2.2.1）；数字位次按 P-25.3.8.1 省略。

本轮只断言稠合母体骨架（P-25.3.2 组装）与取代基位次（P-25.3.3 编号）；
hydro 前缀（P-31.2.2）与指示氢（P-58.2.1）仍由既有通用路径给出，不作断言。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
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


# 近邻负例：已注册稠环与单环烃不得被本轮组分路径改写
NEIGHBORS = [
    ("c1ccc2occc2c1", "1-benzofuran", "1-苯并呋喃"),
    ("C1CCC2CCCCC2C1", "decahydronaphthalene", "十氢萘"),
    ("C1CCc2ccccc2C1", "1,2,3,4-tetrahydronaphthalene", "1,2,3,4-四氢萘"),
    ("C1CCCC1", "cyclopentane", "环戊烷"),
    ("OC1CCCC1", "cyclopentanol", "环戊醇"),
]


@pytest.mark.parametrize("smiles,en,zh", NEIGHBORS)
def test_neighbors_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh
