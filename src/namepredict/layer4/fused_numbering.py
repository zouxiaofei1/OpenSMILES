"""P-25.3.3 稠环编号: 外周骨架编号 + 稠合碳 a/b/c 字母位次 + 准则(a)-(d)。"""
from __future__ import annotations

import re
from collections import Counter, defaultdict

from namepredict.constants import (
    Al, As, B, Bi, Br, C, Cl, F, Ga, Ge, I, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl,
)

# P-25.3.3.1.2(b): F > Cl > Br > I > O > S > Se > Te > N > P > ... > Tl（低位次给更优先杂原子）。
_P145_SENIOR = (F, Cl, Br, I, O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)


def _label_key(lbl) -> tuple[int, str]:
    """locant 串排序键: "4a" → (4, 'a'), "10" → (10, '')。"""
    m = re.match(r"(\d+)([a-z]*)", str(lbl))
    return (int(m.group(1)), m.group(2)) if m else (0, "")


def fused_atoms(rings) -> set[int]:
    """出现在 ≥2 个环的原子(稠合原子)。"""
    counts = defaultdict(int)
    for ring in rings:
        for a in ring:
            counts[a] += 1
    return {a for a, c in counts.items() if c >= 2}


def _top_rings(coords: dict, rings) -> list[int]:
    """环心 (cy 最大, cx 最大) 的最上端环; 平局全保留。"""
    centers = [(r, sum(coords[a][0] for a in ring) / len(ring),
                sum(coords[a][1] for a in ring) / len(ring)) for r, ring in enumerate(rings)]
    cy_max = max(c[2] for c in centers)
    tops = [c for c in centers if abs(c[2] - cy_max) < 1e-9]
    cx_max = max(c[1] for c in tops)
    return [r for r, cx, cy in tops if abs(cx - cx_max) < 1e-9]


def _top_atoms(coords: dict, ring, fused: set[int]) -> list[int]:
    """环内 y 最大(x 平局)的非稠合原子, 平局全保留。"""
    atoms = [a for a in ring if a not in fused]
    if not atoms:
        return []
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
    """沿外部边外边界行走, 用鞋带面积强制顺时针方向。

    稠合系统外边界是由全部外部边组成的单闭环; 每原子恰好 2 个外部边邻居,
    从 start 沿唯一未访问方向走回起点即可。逆时针则整体反转(start 保持首位)。
    """
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
    """P-25.3.3.1.1 全部编号候选 (chain, labels): 起点环/起点原子平局组合。"""
    fused_carbons = {a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C}
    neighbors = _ring_neighbors(rings)
    exterior = _exterior_edges(rings)
    cands: list[tuple[list[int], list[str]]] = []
    for sr in _top_rings(coords, rings):
        starts = _top_atoms(coords, rings[sr], fused)
        if not starts:
            # 最上端环无非稠合原子: 沿顺时针取相邻环中最上端者(P-25.3.3.1.1 兜底)。
            for nb_ring in rings:
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
    locs = sorted((_label_key(labels[chain.index(a)]) for a in atoms if a in chain),
                  key=lambda k: (k[0], k[1]))
    return tuple(locs)


def number_fused_system(mol, rings, coords) -> tuple[list[int], list[str]] | None:
    """P-25.3.3 稠环编号; 依次应用准则(a)-(d) 收窄; 返回 (chain, labels) 或 None。"""
    fused = fused_atoms(rings)
    heteros = _hetero_set(mol, fused) | _hetero_set(mol, set().union(*rings))
    cands = _candidates(mol, rings, coords, fused)
    if not cands:
        return None
    fused_carbons = sorted(a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() == C)
    fused_heteros = sorted(a for a in fused if mol.GetAtomWithIdx(a).GetAtomicNum() != C)
    all_heteros = sorted(heteros)

    def _keep(cands, atoms):
        """保留杂原子/稠合原子位次集合最小的候选。"""
        best = min(_locant_tuples(c[0], c[1], atoms) for c in cands)
        return [c for c in cands if _locant_tuples(c[0], c[1], atoms) == best]

    if all_heteros:
        cands = _keep(cands, all_heteros)  # (a) 低位次给杂原子集合
    # (b) 按 F>Tl 顺序逐元素收窄该元素原子位次
    hetero_by_z = defaultdict(list)
    for a in all_heteros:
        hetero_by_z[mol.GetAtomWithIdx(a).GetAtomicNum()].append(a)
    for z in _P145_SENIOR:
        if z in hetero_by_z:
            cands = _keep(cands, sorted(hetero_by_z[z]))
    if fused_carbons:
        cands = _keep(cands, fused_carbons)  # (c) 低位次给稠合碳
    if fused_heteros:
        cands = _keep(cands, fused_heteros)  # (d) 低位次给稠合杂原子
    return cands[0]
