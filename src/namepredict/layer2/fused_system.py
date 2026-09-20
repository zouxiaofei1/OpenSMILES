"""P-25.3.2.4 稠环拆解为保留母体组分树（fused_info）。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from namepredict.constants import C, P145_SENIOR, P25_SENIOR
from namepredict.layer2.ring_scaffold import (
    _Q,
    component_stem,
    match_fusion_component,
    match_retained,
    omits_fusion_numbers,
    retained_fusion_prefix,
)
from namepredict.layer1.ring_systems import sssr_rings
from namepredict.layer4.numbering_engine import narrow

_TEMPLATE_COUNTS: dict[str, Counter] = {  # 保留模板的元素计数
    sid: Counter(q.GetAtomWithIdx(i).GetAtomicNum() for i in range(q.GetNumAtoms()))
    for sid, q in _Q.items()
}


@dataclass(frozen=True)
class FusedNode:
    """稠环系统的一个组分（母体组分或附加组分，递归）。"""
    scaffold_id: str
    atom_ids: tuple[int, ...]
    ring_indices: frozenset[int]
    fusion_shared: tuple[frozenset[int], ...] = ()
    attached: tuple["FusedNode", ...] = ()
    fused_stem: tuple[str, str] | None = None    # 组分词干 (en, zh)；None = 不可作稠合零件
    fused_prefix: tuple[str, str] | None = None  # 附加组分保留前缀 (en, zh)；None = 走通用规则
    fused_omit_numbers: bool = False             # 稠合描述符省略数字位次（P-25.3.8.1：一级单环烃附加组分）

def _fusion_adj(ring_indices, fusion_edges) -> dict[int, set[int]]:
    """融合图邻接表：环下标 → 相邻环下标（仅两端都在集合内的稠合边）。"""
    adj: dict[int, set[int]] = {r: set() for r in ring_indices}
    for i, j, _ in fusion_edges:
        if i in ring_indices and j in ring_indices:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def _candidates_for(info, rings, fusion_edges, ring_indices) -> dict[frozenset[int], tuple[str, frozenset[int]]]:
    """增长式枚举环集内的保留母体候选（原子集去重、超集剪枝）。"""
    mol = info["mol"]
    adj = _fusion_adj(ring_indices, fusion_edges)
    out: dict[frozenset[int], tuple[str, frozenset[int]]] = {}
    for seed in sorted(ring_indices):
        seed_atoms = frozenset(rings[seed])
        if match_fusion_component(info, seed_atoms) is None:  # P-25.2.1/25.3.2.2.1
            continue
        stack: list[tuple[frozenset[int], frozenset[int]]] = [(seed_atoms, frozenset({seed}))]
        visited: set[frozenset[int]] = set()
        while stack:
            atoms, rset = stack.pop()
            if rset in visited:
                continue
            visited.add(rset)
            sid = match_fusion_component(info, atoms)
            if sid is not None:
                out.setdefault(atoms, (sid, rset))
            for nb in adj:
                if nb in rset or not any(nb in adj[r] for r in rset):
                    continue
                stack.append((atoms | frozenset(rings[nb]), rset | {nb}))
    return out


def _select_base(info, rings, fusion_edges, ring_indices) -> tuple[str, frozenset[int], frozenset[int]] | None:
    """P-25.3.2.4 按准则(a)-(j) 选母体组分，仍并列取环集升序最小。"""
    cands = [(atoms, sid, rset) for atoms, (sid, rset) in
             _candidates_for(info, rings, fusion_edges, ring_indices).items()]
    if not cands:
        return None
    mol = info["mol"]

    def _hetero(c) -> Counter:
        """候选的杂原子计数（排除碳）。"""
        return Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in c[0]
                       if mol.GetAtomWithIdx(i).GetAtomicNum() != C)

    def _key_a(c) -> int:
        """(a) 候选最优先杂原子在 P25_SENIOR 的下标。"""
        h = _hetero(c)
        return min(P25_SENIOR.index(z) for z in h) if h else len(P25_SENIOR)

    cands = narrow(cands, _key_a)  # (a) 含有更优先的杂原子，取键最小
    cands = narrow(cands, lambda c: len(c[2]), reverse=True)   # (b) 环数更多
    cands = narrow(cands, lambda c: tuple(sorted((len(rings[i]) for i in c[2]), reverse=True)),
                   reverse=True)  # (c) 环大小降序最大
    cands = narrow(cands, lambda c: sum(_hetero(c).values()), reverse=True)  # (d) 杂原子总数更多
    cands = narrow(cands, lambda c: len(_hetero(c)), reverse=True)           # (e) 杂原子种类更多

    def _key_f(c):
        """(f) 最高优先杂原子(P145_SENIOR)的 (-rank,计数)。"""
        hc = _hetero(c)
        if not hc:
            return (-len(P145_SENIOR), 0)
        top = min(hc, key=lambda z: P145_SENIOR.index(z))
        return (-P145_SENIOR.index(top), hc[top])

    cands = narrow(cands, _key_f, reverse=True)  # (f) 最高优先性杂原子数更多
    numbering = {c: _numbered_locants(info, rings, fusion_edges, c) for c in cands}  # (g)-(j): 依赖 L4 编号，逐准则收窄
    from namepredict.layer4.fused_numbering import fused_atoms
    from namepredict.layer4.locant_calc import locant_key

    def _locant_tup(c, atoms):
        """取候选下指定原子集的 locant 排序元组。"""
        labels = numbering[c][0]
        return tuple(sorted((locant_key(labels[a]) for a in atoms if a in labels)))

    def _gj(key_fn, *, lowest=False):
        """按 key_fn 对可编号候选取最优值收窄（可编号候选不足 2 个则原样返回）。"""
        scored = [c for c in cands if numbering[c] is not None]
        if len(scored) < 2:
            return cands
        return narrow(scored, key_fn, reverse=not lowest)

    cands = _gj(lambda c: numbering[c][1])  # (g) 水平行环数最多
    cands = _gj(lambda c: _locant_tup(c, [a for a in c[0]
                                          if mol.GetAtomWithIdx(a).GetAtomicNum() != C]), lowest=True)  # (h) 杂原子位次低
    for z in P145_SENIOR:  # (i) 按 F>Tl 逐元素位次低
        cands = _gj(lambda c, z=z: _locant_tup(c, [a for a in c[0]
                                                   if mol.GetAtomWithIdx(a).GetAtomicNum() == z]), lowest=True)

    def _fused_carbons(c):
        """取候选子环集稠合原子的碳原子列表。"""
        fused = fused_atoms([rings[i] for i in sorted(c[2])])
        return [a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C]

    cands = _gj(lambda c: _locant_tup(c, _fused_carbons(c)), lowest=True)  # (j) 稠合碳位次低
    cands = narrow(cands, lambda c: tuple(sorted(c[2]))) 
    atoms, sid, rset = cands[0]
    return sid, atoms, rset


def _numbered_locants(info, rings, fusion_edges, cand):
    """候选母体子环集经 L4 优选取向+P-25.3.3 编号，返回位次表与行数。"""
    _, _, rset = cand
    if len(rset) < 1:
        return None
    sub_rings = [rings[i] for i in sorted(rset)]
    idx_map = {i: k for k, i in enumerate(sorted(rset))}
    sub_edges = [(idx_map[i], idx_map[j], sh) for i, j, sh in fusion_edges
                 if i in rset and j in rset]
    from namepredict.layer4.fused_orientation import preferred_orientations
    from namepredict.layer4.fused_numbering import number_fused_system
    orients = preferred_orientations(info["mol"], sub_rings, sub_edges)
    if not orients:
        return None
    result = number_fused_system(info["mol"], sub_rings, [o.coord_dict() for o in orients])
    if result is None:
        return None
    chain, labels = result
    return dict(zip(chain, labels)), len(orients[0].row)


def _ring_components(remaining, fusion_edges) -> list[frozenset[int]]:
    """剩余环在融合图上的连通分量。"""
    rem = set(remaining)
    adj = _fusion_adj(rem, fusion_edges)
    comps: list[frozenset[int]] = []
    seen: set[int] = set()
    for r in sorted(rem):
        if r in seen:
            continue
        comp: set[int] = set()
        stack = [r]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            comp.add(cur)
            stack.extend(adj[cur] - seen)
        comps.append(frozenset(comp))
    return comps


def _decompose(info, rings, fusion_edges, ring_indices, fusion_shared=()) -> FusedNode | None:
    """递归拆解为 FusedNode 树：选母体后，剩余环按连通分量递归为附加组分。"""
    base = _select_base(info, rings, fusion_edges, ring_indices)
    if base is None:
        return None
    sid, base_atoms, base_rset = base
    attached: list[FusedNode] = []
    for comp in _ring_components(ring_indices - base_rset, fusion_edges):
        comp_shared = tuple(frozenset(s) for i, j, s in fusion_edges
                            if (i in base_rset and j in comp) or (j in base_rset and i in comp))
        node = _decompose(info, rings, fusion_edges, comp, comp_shared)
        if node is not None:
            attached.append(node)
    return FusedNode(
        scaffold_id=sid,
        atom_ids=tuple(sorted(base_atoms)),
        ring_indices=base_rset,
        fusion_shared=tuple(sorted(fusion_shared, key=lambda s: sorted(s))),
        attached=tuple(attached),
        fused_stem=component_stem(sid),
        fused_prefix=retained_fusion_prefix(sid),
        fused_omit_numbers=omits_fusion_numbers(sid),
    )


def decompose_fused_system(info, system) -> FusedNode | None:
    """公共入口: 环系拆解为 FusedNode 树，无候选返回 None。"""
    rings = list(sssr_rings(info["mol"]))
    atom_ids = tuple(system.get("atom_ids") or ())
    if atom_ids:  # 整环系已是一个不可作稠合组分的保留母体（如金刚烷这类笼状桥烃）：稠合拆解无意义且会误判加氢
        sid = match_retained(info, atom_ids)
        if sid is not None and component_stem(sid) is None:
            return None
    node = _decompose(info, rings, system.get("fusion_edges") or [],
                      frozenset(system.get("sssr_indices") or ()))
    return node

