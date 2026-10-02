"""P-23 桥环编号裁决：在 L2 并列候选间按 P-23.3.2 与 P-14.4 选编号。

L2 已给出全部并列最优的 von Baeyer 拆解，此处只做裁决：选中的候选写回
parent["bridged_node"]，返回其位次升序原子表供下游当 chain 用。
"""
from __future__ import annotations

from namepredict.layer4.numbering_engine import (
    _locant_set, _principal_atoms, _unsat_bonds, alpha_locants, by_z,
    narrow, narrow_by_senior, resolve_numbering,
)


def bridged_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回桥环骨架的位次升序原子表；候选不可判定时返回 None。"""
    return resolve_numbering(parent, substituents, "bridged_node",
                             lambda nd: nd.locant_pairs, _narrow_ladder)


def _narrow_ladder(nodes: list, parent: dict, substituents: list, mol, chain: list[int],
                   heteros: list[int]) -> list:
    """依次施加 P-23.3.2.1/.2、P-14.4(c)/(f)/(g)。"""
    key = lambda nd, atoms: _locant_set(nd.numbering, atoms)
    if heteros and mol is not None:  # P-23.3.2.1 集合最低 → .2 逐元素
        # P-23.3.2.2 的序列是 P145_SENIOR 去掉卤素；卤素一价、做不了骨架原子，故等价
        nodes = narrow_by_senior(nodes, key, heteros, by_z(mol, heteros), skip_none=True)
    principal = [a for a in _principal_atoms(parent) if a in chain]
    if principal:  # P-14.4(c) 主特征基团（后缀）位次最低
        nodes = narrow(nodes, lambda nd: key(nd, principal), skip_none=True)
    subs = sorted(s["attach_idx"] for s in (substituents or [])
                  if s.get("attach_idx") in chain)
    if subs:  # P-14.4(f) 取代基位次集合最低
        nodes = narrow(nodes, lambda nd: key(nd, subs), skip_none=True)
    if substituents:  # P-14.4(g) 字母序最前的取代基位次最低
        nodes = narrow(nodes, lambda nd: alpha_locants(nd.numbering, chain, substituents),
                       skip_none=True)
    bonds, _ = _unsat_bonds(parent)
    if bonds:  # P-31.1.4.2 残余平局：复合位次数目最少 → 忽略括号比较 → 全集合最低
        nodes = narrow(nodes, lambda nd: _unsat_key(nd, bonds), skip_none=True)
    return nodes


def _unsat_key(node, bonds) -> tuple | None:
    """多重键的 P-31.1.4.2 键：复合位次数目最少 → 忽略括号内比较 → 全集合最低。"""
    lows, alls, n_comp = [], [], 0
    for a, b in bonds:
        if a not in node.numbering or b not in node.numbering:
            return None
        lo, hi = sorted((node.numbering[a], node.numbering[b]))
        lows.append(lo)
        alls += [lo, hi]
        n_comp += hi - lo != 1
    return (n_comp, tuple(sorted(lows)), tuple(sorted(alls)))
