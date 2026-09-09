"""P-25.3.2.4 稠环拆解: 稠环系统 → 保留母体组分树(fused_info 结构事实)。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from namepredict.constants import (
    Al, As, B, Bi, Br, C, Cl, F, Ga, Ge, I, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl,
)
from namepredict.layer2.ring_scaffold import _Q, match_retained


# P-25.3.2.4(a): N > F > Cl > Br > I > O > S > Se > Te > P > ... > Tl（母体组分选择，N 最优先）。
_P25_SENIOR = (N, F, Cl, Br, I, O, S, Se, Te, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)
# P-25.3.2.4(f): F > Cl > Br > I > O > S > Se > Te > N > P > ... > Tl（同序更高优先杂原子计数，F 最优先）。
_P145_SENIOR = (F, Cl, Br, I, O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)

# 保留模板的元素计数（超集剪枝：当前原子集须 ≤ 某模板计数才可能拼出保留母体）。
_TEMPLATE_COUNTS: dict[str, Counter] = {
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


def _has_template_superset(mol, atom_ids) -> bool:
    """当前原子集元素计数是否 ≤ 某保留模板计数（增长剪枝必要条件）。"""
    current = Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids)
    return any(all(current[z] <= tc.get(z, 0) for z in current) for tc in _TEMPLATE_COUNTS.values())


def _seedable(info, ring_atoms) -> bool:
    """单环是否精确匹配某保留模板（作为增长种子的必要条件）。"""
    return match_retained(info, ring_atoms) is not None


def _candidates_for(info, rings, fusion_edges, ring_indices) -> dict[frozenset[int], tuple[str, frozenset[int]]]:
    """增长式枚举环集内全部保留母体候选 {原子集: (scaffold_id, 环集)}：从可种子环 DFS 并入邻接环，match_retained 精确命中记录、超集剪枝、原子集去重。"""
    mol = info["mol"]
    adj: dict[int, set[int]] = {r: set() for r in ring_indices}
    for i, j, _ in fusion_edges:
        if i in ring_indices and j in ring_indices:
            adj[i].add(j)
            adj[j].add(i)
    out: dict[frozenset[int], tuple[str, frozenset[int]]] = {}
    for seed in sorted(ring_indices):
        seed_atoms = frozenset(rings[seed])
        if not _seedable(info, seed_atoms):
            continue
        stack: list[tuple[frozenset[int], frozenset[int]]] = [(seed_atoms, frozenset({seed}))]
        visited: set[frozenset[int]] = set()
        while stack:
            atoms, rset = stack.pop()
            if rset in visited:
                continue
            visited.add(rset)
            sid = match_retained(info, atoms)
            if sid is not None:
                out.setdefault(atoms, (sid, rset))
            if not _has_template_superset(mol, atoms):
                continue
            for nb in adj:
                if nb in rset or not any(nb in adj[r] for r in rset):
                    continue
                stack.append((atoms | frozenset(rings[nb]), rset | {nb}))
    return out


def _keep_best(cands, key, *, reverse: bool = False):
    """按 key 保留最优候选（max；reverse=True 取 min），每步过滤后候选缩小。"""
    best = (min if reverse else max)((key(c) for c in cands), default=None)
    return [c for c in cands if key(c) == best] if best is not None else cands


def _select_base(info, rings, fusion_edges, ring_indices) -> tuple[str, frozenset[int], frozenset[int]] | None:
    """P-25.3.2.4 按准则(a)-(j) 选母体组分，返回 (scaffold_id, 原子 frozenset, 环 frozenset) 或 None；仍 >1 时取环集升序最小的确定性兜底。"""
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
        """(a) 候选最优先杂原子在 _P25_SENIOR 的下标（无杂原子取最大，N 最优先）。"""
        h = _hetero(c)
        return min(_P25_SENIOR.index(z) for z in h) if h else len(_P25_SENIOR)

    cands = _keep_best(cands, _key_a, reverse=True)  # (a) 含有更优先的杂原子
    cands = _keep_best(cands, lambda c: len(c[2]))   # (b) 环数更多
    cands = _keep_best(cands, lambda c: tuple(sorted((len(rings[i]) for i in c[2]), reverse=True)))  # (c) 环大小降序最大
    cands = _keep_best(cands, lambda c: sum(_hetero(c).values()))  # (d) 杂原子总数更多
    cands = _keep_best(cands, lambda c: len(_hetero(c)))           # (e) 杂原子种类更多

    def _key_f(c):
        """(f) 最高优先杂原子(_P145_SENIOR)的 (-rank, count)；rank 更小优先，同 rank 计数更多优先。"""
        hc = _hetero(c)
        if not hc:
            return (-len(_P145_SENIOR), 0)
        top = min(hc, key=lambda z: _P145_SENIOR.index(z))
        return (-_P145_SENIOR.index(top), hc[top])

    cands = _keep_best(cands, _key_f)  # (f) 最高优先性杂原子数更多
    numbering = {c: _numbered_locants(info, rings, fusion_edges, c) for c in cands}  # (g)-(j): 依赖 L4 优选取向/编号, 逐准则按水平行环数/位次收窄(候选无法编号则跳过该准则)。
    from namepredict.layer4.fused_numbering import fused_atoms
    from namepredict.layer4.locant_key import locant_key

    def _locant_tup(c, atoms):
        """取候选下指定原子集的 locant 排序元组。"""
        labels = numbering[c][0]
        return tuple(sorted((locant_key(labels[a]) for a in atoms if a in labels)))

    def _gj(key_fn, *, reverse=False):
        """按 key_fn 对可编号候选取最优值收窄（候选不足 2 个则原样返回）。"""
        scored = [(c, key_fn(c)) for c in cands if numbering[c] is not None]
        if len(scored) < 2:
            return cands
        best = (max if not reverse else min)(k for _, k in scored)
        return [c for c, k in scored if k == best]

    cands = _gj(lambda c: numbering[c][1])  # (g) 水平行环数最多
    cands = _gj(lambda c: _locant_tup(c, [a for a in c[0]
                                          if mol.GetAtomWithIdx(a).GetAtomicNum() != C]), reverse=True)  # (h) 杂原子位次低
    for z in _P145_SENIOR:  # (i) 按 F>Tl 逐元素位次低
        cands = _gj(lambda c, z=z: _locant_tup(c, [a for a in c[0]
                                                   if mol.GetAtomWithIdx(a).GetAtomicNum() == z]), reverse=True)

    def _fused_carbons(c):
        """取候选子环集稠合原子的碳原子列表。"""
        fused = fused_atoms([rings[i] for i in sorted(c[2])])
        return [a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C]

    cands = _gj(lambda c: _locant_tup(c, _fused_carbons(c)), reverse=True)  # (j) 稠合碳位次低
    cands = _keep_best(cands, lambda c: tuple(sorted(c[2])), reverse=True)  # 兜底：环集升序最小
    atoms, sid, rset = cands[0]
    return sid, atoms, rset


def _numbered_locants(info, rings, fusion_edges, cand):
    """候选母体(子环集)经 L4 优选取向+P-25.3.3 编号，返回 ({原子: locant 串}, 水平行环数) 或 None(该准则跳过)。"""
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
    adj: dict[int, set[int]] = {r: set() for r in rem}
    for i, j, _ in fusion_edges:
        if i in rem and j in rem:
            adj[i].add(j)
            adj[j].add(i)
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
    """递归拆解为 FusedNode 树：选母体组分后，剩余环按连通分量递归为附加组分（fusion_shared 为本组分与父组分的共享原子集）。"""
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
    )


def decompose_fused_system(info, system) -> FusedNode | None:
    """公共入口: 对单个环系拆解为 FusedNode 树（无保留候选返回 None）。"""
    rings = list(info["mol"].GetRingInfo().AtomRings())
    node = _decompose(info, rings, system.get("fusion_edges") or [],
                      frozenset(system.get("sssr_indices") or ()))
    return node

