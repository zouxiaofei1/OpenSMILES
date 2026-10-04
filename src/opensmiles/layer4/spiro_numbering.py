"""P-24 螺环编号：单环组分螺描述符候选 + 组分式螺环逐组分编号。
单环组分按 P-24.2.2/2.3、P-24.2.4.1.2 收窄候选。
组分式螺环按 P-24.5.2 让螺稠合位次最小、第 k 组分带 k 撇号。
"""
from __future__ import annotations

from dataclasses import replace
from itertools import product

from opensmiles.constants import C
from opensmiles.layer4.locant_calc import locant_key
from opensmiles.layer4.numbering_engine import (
    _bond_locant_pairs, _locant_set, _principal_atoms, _unsat_bonds, alpha_locants, by_z,
    candidates, hetero_atoms, narrow, narrow_by_senior, resolve_numbering,
)

APOSTROPHE = "'"
_MAX_COMBOS = 4096  # 组分候选组合上限


def spiro_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回螺环骨架的位次升序原子表；候选不可判定时返回 None。"""
    return resolve_numbering(parent, substituents, "spiro_node",
                             lambda nd: nd.descriptor_superscripts, _narrow_ladder)


def _narrow_ladder(nodes: list, parent: dict, substituents: list, mol, chain: list[int],
                   heteros: list[int]) -> list:
    """依次施加 P-24.2.2/2.3、P-24.2.4.1.2、P-14.4。"""
    key = lambda nd, atoms: _locant_set(nd.numbering, atoms)
    nodes = narrow(nodes, lambda nd: key(nd, list(nd.free_spiro_atoms)))  # P-24.2.2.1 / P-24.2.3.1 螺原子位次集合最低
    nodes = narrow(nodes, lambda nd: (nd.descriptor, nd.descriptor_superscripts))  # P-24.2.2.2 / P-24.2.3.2 描述符数字取小
    if heteros and mol is not None:  # P-24.2.4.1.2(a) 集合 → (b) 逐元素
        nodes = narrow_by_senior(nodes, key, heteros, by_z(mol, heteros), skip_none=True)
    principal = [a for a in _principal_atoms(parent) if a in chain]
    if principal:  # P-14.4(c) 主特征基团（后缀）位次最低
        nodes = narrow(nodes, lambda nd: key(nd, principal), skip_none=True)
    bonds, doubles = _unsat_bonds(parent)
    if bonds and mol is not None:  # P-14.4(e) 双键位次最低
        nodes = narrow(nodes, lambda nd: _bond_locant_pairs(nd.numbering, bonds, doubles), skip_none=True)
    subs = sorted(s["attach_idx"] for s in (substituents or [])
                  if s["attach_idx"] in chain)
    if subs:  # P-14.4(f) 取代基位次集合最低
        nodes = narrow(nodes, lambda nd: key(nd, subs), skip_none=True)
    if substituents:  # P-14.4(g) 字母序最前的取代基位次最低
        nodes = narrow(nodes, lambda nd: alpha_locants(nd.numbering, chain, substituents),
                       skip_none=True)
    return nodes


def fbs_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回组分式螺环骨架的位次升序原子表；候选不可判定时返回 None。"""
    nodes = candidates(parent, "fbs_node")
    mol = parent.get("mol")
    if not nodes or mol is None:
        return None
    node = nodes[0]
    cand_lists = [_shrink(comp, mol, parent, substituents) for comp in node.components]
    if any(not c for c in cand_lists):
        return None
    chosen = _joint_pick(node, cand_lists)
    if chosen is None:
        return None
    chain: list[int] = []
    labels: list[str] = []
    seen: set[int] = set()
    for comp, num in zip(node.components, chosen):
        marks = APOSTROPHE * comp.index  # P-24.6：第 k 个组分位次带 k 个撇号
        for atom, lbl in zip(num.chain, num.labels):
            if atom in seen:
                continue  # 螺原子同属两组分：全局链每原子只留一次（取先列组分的位次）
            seen.add(atom)
            chain.append(atom)
            labels.append(str(lbl) + marks)
    comps = tuple(replace(c, numberings=(num,)) for c, num in zip(node.components, chosen))
    parent["fbs_node"] = replace(node, components=comps)  # L5 依它取组分名与连接点
    parent["numbering_scaffold"] = {"scaffold_id": "fused_bridged_spiro", "labels": tuple(labels)}
    return chain


def _fbs_locant_set(num, atoms) -> tuple:
    """指定原子集在组分候选下的位次集合。"""
    loc = num.locants
    return tuple(sorted(locant_key(loc[a]) for a in atoms if a in loc))


def _fbs_feature_key(mol, comp, substituents):
    """组分候选的命名特征：杂原子/螺原子/取代基位次（对称等价者相同）。"""
    def feat(num) -> tuple:
        loc = num.locants
        het = tuple(sorted((mol.GetAtomWithIdx(a).GetAtomicNum(), str(loc[a]))
                           for a in comp.atom_ids if a in loc
                           and mol.GetAtomWithIdx(a).GetAtomicNum() != C))
        subs = tuple(sorted(str(loc[s["attach_idx"]]) for s in (substituents or ())
                            if s.get("attach_idx") in loc))
        return (het, _fbs_locant_set(num, comp.spiro_atoms), subs)
    return feat


def _shrink(comp, mol, parent, substituents) -> list:
    """单组分内收窄编号候选：P-14.4(c) 后缀→螺位次→杂原子→取代基。"""
    cands = list(comp.numberings)
    if not cands:
        return []
    principal = [a for a in _principal_atoms(parent) if a in comp.atom_ids]
    if principal:  # P-14.4(c)：主特征基团（后缀）位次最低
        cands = narrow(cands, lambda c: _fbs_locant_set(c, principal), skip_none=True)
    cands = narrow(cands, lambda c: _fbs_locant_set(c, comp.spiro_atoms), skip_none=True)  # P-24.5.2/4：螺稠合位次优先于 'a' 前缀
    heteros = hetero_atoms(mol, comp.atom_ids)
    if heteros:  # P-14.4(e)/P-22.2.2.1.3：杂原子位次集合最低
        cands = narrow_by_senior(cands, _fbs_locant_set, heteros, by_z(mol, heteros),
                                 skip_none=True)
    subs = sorted({s["attach_idx"] for s in (substituents or ())
                   if s.get("attach_idx") in comp.atom_ids})
    if subs:  # P-14.4(f)：取代基位次集合最低
        cands = narrow(cands, lambda c: _fbs_locant_set(c, subs), skip_none=True)
    return cands


def _listing_key(node, combo) -> tuple | None:
    """按名称列出顺序展平的螺位次元组（P-24.6.1：列出顺序低者在前）。"""
    out = []
    for s, a, b in node.links:
        la, lb = combo[a].locants.get(s), combo[b].locants.get(s)
        if la is None or lb is None:
            return None
        out.append(locant_key(la))
        out.append(locant_key(lb))
    return tuple(out)


def _joint_pick(node, cand_lists: list):
    """引用序位次元组最低者胜出；并列取先（P-24.6.1）。"""
    total = 1
    for c in cand_lists:
        total *= len(c)
    if total > _MAX_COMBOS:
        return None
    best, best_key = None, None
    for combo in product(*cand_lists):
        key = _listing_key(node, combo)
        if key is None:
            continue
        if best_key is None or key < best_key:
            best, best_key = combo, key
    return best
