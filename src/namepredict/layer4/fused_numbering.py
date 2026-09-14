"""P-25.3.3 稠环编号：外周序、稠合位次、(a)-(f)、CIP 破局。"""
from __future__ import annotations

from collections import Counter, defaultdict

from namepredict.constants import C
from namepredict.layer4.indicated_hydrogen import saturated_ring_atoms
from namepredict.layer4.locant_calc import locant_key
from namepredict.layer4.numbering_engine import narrow, narrow_by_senior

INDICATED_H = object()  # sub_layers 哨兵层：P-25.3.3.1.2(f) 指示氢最小化。


def fused_atoms(rings) -> set[int]:
    """出现在 ≥2 个环的原子(稠合原子)。"""
    counts = defaultdict(int)
    for ring in rings:
        for a in ring:
            counts[a] += 1
    return {a for a, c in counts.items() if c >= 2}


def _top_rings(coords: dict, rings) -> list[int]:
    """最上端水平行最右环(P-25.3.3.1.1): y 聚行后取 x 最大。"""
    centers = [(r, sum(coords[a][0] for a in ring) / len(ring),
                sum(coords[a][1] for a in ring) / len(ring)) for r, ring in enumerate(rings)]
    span = max(max(coords[a][1] for a in ring) - min(coords[a][1] for a in ring) for ring in rings)
    cy_max = max(c[2] for c in centers)
    tops = [c for c in centers if cy_max - c[2] <= 0.25 * span]  # 视为同处一个水平行
    cx_max = max(c[1] for c in tops)
    return [r for r, cx, cy in tops if abs(cx - cx_max) < 1e-9]


def _top_atoms(coords: dict, ring, fused: set[int]) -> list[int]:
    """环内 y 最大的非稠合起点原子(P-25.3.3.1.1), 平局全保留。"""
    atoms = [a for a in ring if a not in fused]
    if not atoms:
        return []
    neighbors = _ring_neighbors([ring])
    adjacent = [a for a in atoms if any(nb in fused for nb in neighbors[a])]
    if adjacent:  # 起点取偏环顶端的顶点会整体错位一位。
        atoms = adjacent
    top = max(atoms, key=lambda a: (coords[a][1], coords[a][0]))
    return [a for a in atoms
            if abs(coords[a][1] - coords[top][1]) < 1e-9 and abs(coords[a][0] - coords[top][0]) < 1e-9]


def _ring_neighbors(rings) -> dict[int, set[int]]:
    """原子 → 全部环内键邻居。"""
    neigh: dict[int, set[int]] = defaultdict(set)
    for ring in rings:
        n = len(ring)
        for k, a in enumerate(ring):
            neigh[a].add(ring[(k + 1) % n])
            neigh[a].add(ring[(k - 1) % n])
    return {a: ns for a, ns in neigh.items()}


def _exterior_edges(rings) -> frozenset[frozenset[int]]:
    """外部边: 只属于一个环的边(非稠合共享边); 稠合系统外边界由其组成。"""
    counts: Counter = Counter()
    for ring in rings:
        n = len(ring)
        for k, a in enumerate(ring):
            counts[frozenset((a, ring[(k + 1) % n]))] += 1
    return frozenset(e for e, c in counts.items() if c == 1)


def _boundary_walk(coords: dict, neighbors: dict, exterior: frozenset[frozenset[int]],
                   start: int) -> list[int]:
    """沿外部边单闭环行走外边界，鞋带面积强制顺时针；逆时针则整体反转。"""
    walk = [start]
    prev, cur = None, start
    while True:
        ext = [nb for nb in neighbors[cur]
               if frozenset((cur, nb)) in exterior and nb != prev]
        if not ext:
            break
        nxt = ext[0]
        if nxt == start:
            break
        walk.append(nxt)
        prev, cur = cur, nxt
    area = 0.0
    for i, a in enumerate(walk):
        x1, y1 = coords[a]
        x2, y2 = coords[walk[(i + 1) % len(walk)]]
        area += x1 * y2 - x2 * y1
    if area > 0:  # 数学正方向(逆时针) → 反转成顺时针
        walk = [walk[0]] + list(reversed(walk[1:]))
    return walk


def _assign_labels(walk: list[int], fused_carbons: set[int], mol) -> tuple[list[int], list[str]]:
    """非稠合原子/稠合杂原子→下一数字; 稠合碳→紧邻前数字+a/b/c 递增。"""
    chain: list[int] = []
    labels: list[str] = []
    num, letter, last_num = 0, 0, 0
    for atom in walk:
        if atom in fused_carbons:
            letter += 1
            chain.append(atom)
            labels.append(f"{last_num}{chr(96 + letter)}")
        else:
            num += 1
            letter = 0
            last_num = num
            chain.append(atom)
            labels.append(str(num))
    return chain, labels


def _hetero_set(mol, atoms: set[int]) -> set[int]:
    """原子集中非碳原子。"""
    return {a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != C}


def _candidates(mol, rings, coords, fused: set[int]) -> list[tuple[list[int], list[str]]]:
    """P-25.3.3.1.1 全部编号候选: 起点环与起点原子平局组合。"""
    fused_carbons = {a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C}
    neighbors = _ring_neighbors(rings)
    exterior = _exterior_edges(rings)
    cands: list[tuple[list[int], list[str]]] = []
    for sr in _top_rings(coords, rings):
        starts = _top_atoms(coords, rings[sr], fused)
        if not starts:
            for nb_ring in rings:  # 最上端环无起点时取相邻环者(P-25.3.3.1.1 兜底)
                if set(nb_ring) & set(rings[sr]):
                    starts = _top_atoms(coords, nb_ring, fused)
                    if starts:
                        break
        for s in starts:
            walk = _boundary_walk(coords, neighbors, exterior, s)
            chain, labels = _assign_labels(walk, fused_carbons, mol)
            if chain:
                cands.append((chain, labels))
    return cands


def _locant_tuples(chain: list[int], labels: list[str], atoms: list[int]) -> tuple:
    """候选编号下某原子集的位次元组(按 locant 键排序)。"""
    locs = sorted((locant_key(labels[chain.index(a)]) for a in atoms if a in chain),
                  key=lambda k: (k[0], k[1]))
    return tuple(locs)


def number_fused_system(mol, rings, coords: list[dict], sub_layers=None,
                        alpha_subs=None) -> tuple[list[int], list[str]] | None:
    """P-25.3.3 稠环编号: 依(a)-(d) 收窄，返回候选或 None。"""
    fused = fused_atoms(rings)
    heteros = _hetero_set(mol, fused) | _hetero_set(mol, set().union(*rings))
    cands: list[tuple[list[int], list[str]]] = []
    for c in coords:
        cands.extend(_candidates(mol, rings, c, fused))
    if not cands:
        return None
    fused_carbons = sorted(a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C)
    fused_heteros = sorted(a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() != C)
    all_heteros = sorted(heteros)

    if all_heteros:  # (a) 低位次给杂原子集合 (b) 按 F>Tl 顺序逐元素收窄该元素原子位次
        hetero_by_z = defaultdict(list)
        for a in all_heteros:
            hetero_by_z[mol.GetAtomWithIdx(a).GetAtomicNum()].append(a)
        cands = narrow_by_senior(cands, lambda c, at: _locant_tuples(c[0], c[1], at),
                                 all_heteros, hetero_by_z)
    if fused_carbons:
        cands = narrow(cands, lambda c: _locant_tuples(c[0], c[1], fused_carbons))  # (c) 低位次给稠合碳
    if fused_heteros:
        cands = narrow(cands, lambda c: _locant_tuples(c[0], c[1], fused_heteros))  # (d) 低位次给稠合杂原子
    ring_atoms = set().union(*rings)
    ind_h_sats = sorted(saturated_ring_atoms(mol, ring_atoms))  # (f) 指示氢候选位：环内仅以单键连邻环原子且带 H 的饱和位

    def _as_indicated(cands):
        """P-25.3.3.1.2(f)：指示氢位次最小化；偶数个饱和带氢位时不适用。"""
        if not len(ind_h_sats) % 2:
            return cands

        def _key(c):
            """候选的指示氢位次键：全部带氢环位的位次集合升序；位次不全在链内返回 None。"""
            locs = _locant_tuples(c[0], c[1], ind_h_sats)
            return locs if len(locs) == len(ind_h_sats) else None
        return narrow(cands, _key, skip_none=True)

    for layer in (sub_layers or ()):
        if len(cands) <= 1:
            break
        if layer is INDICATED_H:
            cands = _as_indicated(cands)  # (f) 指示氢位次
        elif layer:
            cands = narrow(cands, lambda c: _locant_tuples(c[0], c[1], sorted(layer)))  # 镜像平局: 逐层按位次集合最小化收窄
    if len(cands) > 1 and alpha_subs:  # P-14.5: 位次集合仍相同时，字母序最前的取代基得最低位次
        def _alpha_key(c):
            """字母序键：[(取代基字母序键, 其位次键)] 排序元组。"""
            return tuple(sorted((k, locant_key(c[1][c[0].index(a)]))
                                for k, a in alpha_subs if a in c[0]))
        cands = narrow(cands, _alpha_key)
    if len(cands) > 1:  # P-14.4(j)：位次准则全平局时按 CIP 描述符定方向。
        from namepredict.layer4.numbering_engine import _chain_rs_codes, _rs_locant_key
        codes = _chain_rs_codes(mol, cands[0][0])
        if codes:
            cands = narrow(cands, lambda c: _rs_locant_key(codes, c[0], c[1]))
    return cands[0]
