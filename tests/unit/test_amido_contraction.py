# IUPAC: P-66.1.1.4.3
# Layer: L3,L5
"""R-C(=O)-NH-* 残基的 amido 收缩（P-66.1.1.4.3 方法 1, PIN）：
acetyl/formyl/benzoyl 三词收成 acetamido/formamido/benzamido（gold/ChEBI 全量仅此三词；
长链/烯酰/被取代苯甲酰/杂环羰酰 gold 保持 acylamino 方法 2——见 test_acyl_oyl.py、
test_ring_acyl_carbonyl.py 锁定，本规则不得误伤）。中文 gold 对糖类用「乙酰氨基」与
本表不完全对齐，故仅精确断言英文。

生成点：L5 assembler._mononuclear_radical_names 单取代 azane 分支（现 free_to_yl 出
acetylamino）；括号由 as_substituent._radical_yl_from_sub composite 判定给 (…amino)。
"""

from __future__ import annotations

import pytest

from namepredict.constants import normalize_en
from namepredict.namer import SMILESNNamer

# 片段级正例：* 锚在 N 的 -NH-C(=O)-R 残基，单 N-酰基且 N 有 H。
POS_FRAG = [
    ("*NC(=O)C", "acetamido"),
    ("*NC=O", "formamido"),
    ("*NC(=O)c1ccccc1", "benzamido"),
]


@pytest.mark.parametrize("smiles,en", POS_FRAG)
def test_amido_fragment(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 片段级近邻负例：不在收缩表内不得误伤（长链酰方法 2、被取代苯甲酰仍 acylamino）。
NEG_FRAG = [
    ("*NC(=O)CCC", "butanoylamino"),
    ("*NC(=O)c1ccc(Cl)cc1", "4-chlorobenzoylamino"),
    ("*NC(=O)c1ccco1", "furan-2-carbonylamino"),
]


@pytest.mark.parametrize("smiles,en", NEG_FRAG)
def test_amido_fragment_negative(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例（含 N,N 双乙酰…等倍增场景落 di- 而非 bis-）。
POS_WHOLE = [
    ("CC(O)=N[C@@H](Cc1ccccc1)C(=O)O", "(2S)-2-acetamido-3-phenylpropanoic acid"),
    ("CC(O)=Nc1ccc(C(=O)O)c(O)c1", "4-acetamido-2-hydroxybenzoic acid"),
    ("O=C(O)CN=CO", "2-formamidoacetic acid"),
    ("O=C(O)CN=C(O)c1ccccc1", "2-benzamidoacetic acid"),
    ("CC(O)=Nc1c(I)c(N=C(C)O)c(I)c(C(=O)O)c1I",
     "3,5-diacetamido-2,4,6-triiodobenzoic acid"),
]


@pytest.mark.parametrize("smiles,en", POS_WHOLE)
def test_amido_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
