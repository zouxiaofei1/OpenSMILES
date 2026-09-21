"""P-24.5~24.7 组分式螺环编号：逐组分定编号与全局带撇位次表。

组分自身编号由 L2 给的候选（保留名固定编号 / P-25.3.3 / P-23）在 L4 收窄：
先按 P-24.5.2 / P-24.5.4 让螺稠合位次取最小，再要求并列候选的命名特征全同。
第 k 个组分的位次统一带 k 个撇号，写进 numbering_scaffold.labels，供取代基、
后缀与立体位次沿用。
"""
from __future__ import annotations

from dataclasses import replace
from itertools import product

from namepredict.layer4.locant_calc import locant_key
from namepredict.layer4.numbering_engine import _principal_atoms, narrow, narrow_by_senior

CARBON = 6
APOSTROPHE = "'"
_MAX_COMBOS = 4096  # 组分候选组合上限


def fbs_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回组分式螺环骨架的位次升序原子表；候选不可判定时返回 None。"""
    nodes = _candidates(parent)
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


def _candidates(parent: dict) -> list:
    """L2 下传的并列节点；无组分式节点返回空表。"""
    nodes = parent.get("fbs_nodes")
    if nodes:
        return list(nodes)
    node = parent.get("fbs_node")
    return [node] if node is not None else []


def _spiro_set(num, spiros) -> tuple:
    """组分内螺原子本位次集合（按位次排序）。"""
    loc = num.locants
    return tuple(sorted(locant_key(loc[a]) for a in spiros if a in loc))


def _locant_set(num, atoms) -> tuple:
    """指定原子集在组分候选下的位次集合。"""
    loc = num.locants
    return tuple(sorted(locant_key(loc[a]) for a in atoms if a in loc))


def _hetero_atoms(mol, comp) -> list[int]:
    """组分的骨架杂原子（按原子序表序）。"""
    return [a for a in comp.atom_ids if mol.GetAtomWithIdx(a).GetAtomicNum() != CARBON]


def _feature_key(mol, comp, substituents):
    """组分候选的命名特征：杂原子/螺原子/取代基位次（对称等价者相同）。"""
    def feat(num) -> tuple:
        loc = num.locants
        het = tuple(sorted((mol.GetAtomWithIdx(a).GetAtomicNum(), str(loc[a]))
                           for a in comp.atom_ids if a in loc
                           and mol.GetAtomWithIdx(a).GetAtomicNum() != CARBON))
        subs = tuple(sorted(str(loc[s["attach_idx"]]) for s in (substituents or ())
                            if s.get("attach_idx") in loc))
        return (het, _spiro_set(num, comp.spiro_atoms), subs)
    return feat


def _shrink(comp, mol, parent, substituents) -> list:
    """单个组分内收窄编号候选：P-14.4(c) 后缀 → 螺位次 → 杂原子 → 取代基。

    对称保留母体的自同构（如 2-benzofuran 的 1/3 互换）在这里无法分辨，
    剩下的并列交给 _joint_pick 按引用序列位次元组裁决（P-24.6.1）。
    """
    cands = list(comp.numberings)
    if not cands:
        return []
    principal = [a for a in _principal_atoms(parent) if a in comp.atom_ids]
    if len(cands) > 1 and principal:  # P-14.4(c)：主特征基团（后缀）位次最低
        cands = narrow(cands, lambda c: _locant_set(c, principal), skip_none=True)
    if len(cands) > 1:  # P-24.5.2 / P-24.5.4：螺稠合位次优先于 'a' 前缀位次
        cands = narrow(cands, lambda c: _spiro_set(c, comp.spiro_atoms), skip_none=True)
    heteros = _hetero_atoms(mol, comp)
    if len(cands) > 1 and heteros:  # P-14.4(e)/P-22.2.2.1.3：杂原子位次集合最低
        by_z: dict[int, list[int]] = {}
        for a in heteros:
            by_z.setdefault(mol.GetAtomWithIdx(a).GetAtomicNum(), []).append(a)
        cands = narrow_by_senior(cands, _locant_set, heteros, by_z, skip_none=True)
    subs = sorted({s["attach_idx"] for s in (substituents or ())
                   if s.get("attach_idx") in comp.atom_ids})
    if len(cands) > 1 and subs:  # P-14.4(f)：取代基位次集合最低
        cands = narrow(cands, lambda c: _locant_set(c, subs), skip_none=True)
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
    """引用序位次元组最低的组合胜出；并列取先（P-24.6.1 / P-24.7.2）。"""
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

