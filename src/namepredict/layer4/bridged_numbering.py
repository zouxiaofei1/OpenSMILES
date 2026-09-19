"""P-23 桥环编号裁决：在 L2 并列候选间按 P-23.3.2 与 P-14.4 选编号。

L2 已给出全部并列最优的 von Baeyer 拆解，此处只做裁决：选中的候选写回
parent["bridged_node"]，返回其位次升序原子表供下游当 chain 用。
"""
from __future__ import annotations

from namepredict.layer4.numbering_engine import _locant_set, _principal_atoms, narrow
from namepredict.layer4.numbering_engine import narrow_by_senior
from namepredict.tools.re import alpha_order_key

CARBON = 6


def bridged_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回桥环骨架的位次升序原子表；候选不可判定时返回 None。"""
    nodes = _candidates(parent)
    chain = list(parent.get("chain") or ())
    if not nodes or not chain:
        return None
    mol = parent.get("mol")
    heteros = _hetero_atoms(mol, chain)
    nodes = _narrow_ladder(nodes, parent, substituents, mol, chain, heteros)
    best = _pick_equivalent(nodes, _feature_key(mol, chain, parent, substituents))
    if best is None:
        return None
    parent["bridged_node"] = best  # L5 依它取描述符与上标，须与选中的编号自洽
    return sorted(best.numbering, key=best.numbering.get)


def _candidates(parent: dict) -> list:
    """L2 下传的并列候选；无桥环节点返回空表。"""
    nodes = parent.get("bridged_nodes")
    if nodes:
        return list(nodes)
    node = parent.get("bridged_node")
    return [node] if node is not None else []


def _hetero_atoms(mol, chain: list[int]) -> list[int]:
    """骨架上的杂原子，按 chain 序。"""
    if mol is None:
        return []
    return [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != CARBON]


def _by_z(mol, atoms: list[int]) -> dict[int, list[int]]:
    """原子序数 → 该元素的骨架原子列表。"""
    out: dict[int, list[int]] = {}
    for a in atoms:
        out.setdefault(mol.GetAtomWithIdx(a).GetAtomicNum(), []).append(a)
    return out


def _narrow_ladder(nodes: list, parent: dict, substituents: list, mol, chain: list[int],
                   heteros: list[int]) -> list:
    """依次施加 P-23.3.2.1/.2、P-14.4(c)/(f)/(g)。"""
    key = lambda nd, atoms: _locant_set(nd.numbering, atoms)
    if heteros and mol is not None:  # P-23.3.2.1 集合最低 → .2 逐元素
        # P-23.3.2.2 的序列是 P145_SENIOR 去掉卤素；卤素一价、做不了骨架原子，故等价
        nodes = narrow_by_senior(nodes, key, heteros, _by_z(mol, heteros), skip_none=True)
    principal = [a for a in _principal_atoms(parent) if a in chain]
    if principal:  # P-14.4(c) 主特征基团（后缀）位次最低
        nodes = narrow(nodes, lambda nd: key(nd, principal), skip_none=True)
    subs = sorted(s["attach_idx"] for s in (substituents or [])
                  if s.get("attach_idx") in chain)
    if subs:  # P-14.4(f) 取代基位次集合最低
        nodes = narrow(nodes, lambda nd: key(nd, subs), skip_none=True)
    if len(nodes) > 1 and substituents:  # P-14.4(g) 字母序最前的取代基位次最低
        nodes = narrow(nodes, lambda nd: _alpha_locants(nd, chain, substituents), skip_none=True)
    return nodes


def _alpha_locants(node, chain: list[int], substituents: list) -> tuple:
    """按前修饰基字母序排列的位次元组（P-14.4(g)）。"""
    pairs = sorted(
        (alpha_order_key(s.get("en") or ""), node.numbering.get(s["attach_idx"], 0))
        for s in substituents if s.get("attach_idx") in chain)
    return tuple(loc for _, loc in pairs)


def _feature_key(mol, chain: list[int], parent: dict, substituents: list):
    """候选渲染特征：描述符 + 上标 + 杂原子/后缀/取代基位次。"""
    def feat(node) -> tuple:
        het = tuple(sorted((mol.GetAtomWithIdx(a).GetAtomicNum(), node.numbering[a])
                           for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != CARBON)) \
            if mol is not None else ()
        suffix = tuple(sorted(node.numbering[a] for a in _principal_atoms(parent) if a in chain))
        subs = tuple(sorted(node.numbering[s["attach_idx"]] for s in (substituents or [])
                            if s.get("attach_idx") in chain))
        return (node.descriptor, node.locant_pairs, het, suffix, subs)
    return feat


def _pick_equivalent(nodes: list, feat):
    """并列候选特征全同才取首个，否则视为不可判定。"""
    if not nodes:
        return None
    return nodes[0] if len({feat(nd) for nd in nodes}) == 1 else None
