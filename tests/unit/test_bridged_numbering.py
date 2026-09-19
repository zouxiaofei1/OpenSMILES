# IUPAC: P-23.3.2 / P-14.4
# Layer: L4
#
# 桥环编号裁决：L2 留下并列的 BridgedNode 候选，L4 按 P-23.3.2（杂原子位次）
# 与 P-14.4(c)(f)(g)（后缀/取代基位次、字母序）选出唯一编号。
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_select import select_parent
from namepredict.layer4.numbering_engine import _fused_numbering, orient_numbering

# 每条 = (SMILES, 期望的杂原子符号→位次)。均为并列候选的实际分歧点。
HETERO_CASES = [
    # 8 个并列候选 N 都在 3 位：位次集合相同，不得改变结果
    ("C12CNCC(CC1)CC2", {"N": 3}),
    # 候选分歧为 N=1/N=4：P-23.3.2.1 取最小集合，N 应得 1
    ("OC1CN2CCC1CC2", {"N": 1}),
    # 候选给出 {S=2,N=5} 与 {S=4,N=1}：最小集合 {1,4} 胜（金标 4-thia-1-azabicyclo[3.2.0]hept-2-en-7-one）
    ("O=C1C[C@H]2SC=CN12", {"S": 4, "N": 1}),
]


def _bridged_parent(smiles: str) -> dict:
    """取该分子唯一的桥环 parent。"""
    mol = preprocess(smiles)
    assert mol is not None
    parents = [p for p in select_parent(analyze(mol)) if p.get("scaffold_id") == "bridged"]
    assert len(parents) == 1, smiles
    return parents[0]


@pytest.mark.parametrize("smiles,want", HETERO_CASES)
def test_heteroatom_locants(smiles, want):
    """P-23.3.2.1：杂原子位次集合最低。"""
    mol = preprocess(smiles)
    chain = orient_numbering(_bridged_parent(smiles), [])
    assert chain is not None
    got = {mol.GetAtomWithIdx(a).GetSymbol(): chain.index(a) + 1
           for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6}
    assert got == want


@pytest.mark.parametrize("smiles", [c[0] for c in HETERO_CASES])
def test_chain_is_locant_order_and_covers_skeleton(smiles):
    """L4 返回的 chain 按位次升序，且恰好覆盖骨架原子一次。"""
    parent = _bridged_parent(smiles)
    chain = orient_numbering(parent, [])
    assert chain is not None
    assert sorted(chain) == sorted(parent["chain"])
    assert len(set(chain)) == len(chain)


def test_suffix_locant_decides_among_hetero_ties():
    """P-14.4(c)：杂原子位次打平后由主特征基团（后缀）位次决定。"""
    smiles = "OC1CN2CCC1CC2"  # 金标 1-azabicyclo[2.2.2]octan-3-ol
    mol = preprocess(smiles)
    o = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    oh_c = next(n.GetIdx() for n in mol.GetAtomWithIdx(o).GetNeighbors())  # 后缀锚点是带 OH 的碳
    parent = _bridged_parent(smiles)
    chain = orient_numbering(parent, [{"attach_idx": oh_c, "en": "hydroxy", "zh": "羟基"}])
    assert chain is not None
    assert chain.index(oh_c) + 1 == 3


@pytest.mark.parametrize("smiles", ["C1C2CC3CC1CC(C2)C3", "C12CNCC(CC1)CC2"])
def test_numbering_is_deterministic(smiles):
    """并列候选渲染特征全同时反复调用结果一致，不随候选顺序漂移。"""
    parent = _bridged_parent(smiles)
    assert orient_numbering(parent, []) == orient_numbering(parent, [])


def test_fused_numbering_refuses_bridged():
    """桥环不得走 P-25.3.3 稠环外周编号（chain_fused 判据对桥环恒真）。"""
    parent = _bridged_parent("C1CC2CCC1C2")
    assert _fused_numbering(parent, list(parent["chain"])) is None


def test_non_bridged_parent_unaffected():
    """非桥环 parent 不得被桥环分支接管。"""
    for smiles in ("c1ccccc1", "C1CCCCC1", "c1ccc2ccccc2c1"):
        mol = preprocess(smiles)
        parent = select_parent(analyze(mol))[0]
        assert parent.get("bridged_node") is None, smiles
        chain = orient_numbering(parent, [])
        assert chain is not None and sorted(chain) == sorted(parent["chain"])
