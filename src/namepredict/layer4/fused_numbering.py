"""P-25.3.3 稠环编号: 外周骨架编号 + 稠合碳 a/b/c 字母位次 + 准则(a)-(d)。"""
from __future__ import annotations

import re
from collections import Counter, defaultdict

from namepredict.constants import (
    Al, As, B, Bi, Br, C, Cl, F, Ga, Ge, I, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl,
)

# P-25.3.3.1.2(b): F > Cl > Br > I > O > S > Se > Te > N > P > ... > Tl（低位次给更优先杂原子）。
_P145_SENIOR = (F, Cl, Br, I, O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)

# 双 N 对称 1,3-/1,2-二唑: N1 须吡咯型 N(价 3 带取代/带 H), 排除双 N 互换的镜像匹配。
_DIAZOLE_IDS = ("imidazole", "pyrazole")


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
    """最上端"水平行"内最右环(P-25.3.3.1.1): 同一水平行内各环同属最上端, 平局取最右。6-5-6 行的变形五元环环心被顶点抬偏 ±0.18, 单比环心 y 会把中间环误判成最上端环(起点落到中间环唯一非稠合原子上), 故先按环心 y 聚成行(阈值为最大环高的 1/4)再取行内 x 最大者; 平局全保留。"""
    centers = [(r, sum(coords[a][0] for a in ring) / len(ring),
                sum(coords[a][1] for a in ring) / len(ring)) for r, ring in enumerate(rings)]
    span = max(max(coords[a][1] for a in ring) - min(coords[a][1] for a in ring) for ring in rings)
    cy_max = max(c[2] for c in centers)
    tops = [c for c in centers if cy_max - c[2] <= 0.25 * span]  # 视为同处一个水平行
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
    """沿外部边单闭环行走外边界（每原子恰 2 个外部边邻居，从 start 沿唯一未访问方向走回起点），鞋带面积强制顺时针；逆时针整体反转、start 保持首位。"""
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
            for nb_ring in rings:  # 最上端环无非稠合原子: 沿顺时针取相邻环中最上端者(P-25.3.3.1.1 兜底)
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


def number_fused_system(mol, rings, coords, sub_layers=None,
                        alpha_subs=None) -> tuple[list[int], list[str]] | None:
    """P-25.3.3 稠环编号: 依准则(a)-(d) 收窄; 返回 (chain, labels) 或 None。coords 可为单 dict 或平局候选列表，一并枚举跨镜像收窄。sub_layers 为按优先级排列的环外附着原子组（纯碳环上 (a)-(d) 全平局，须逐层做位次集合最小化收窄，否则编号方向随候选枚举顺序漂移）；alpha_subs 为 [(字母序键, 附着原子)]，位次集合仍相同时按 P-14.5 把最低位次给字母序最前者。"""
    fused = fused_atoms(rings)
    heteros = _hetero_set(mol, fused) | _hetero_set(mol, set().union(*rings))
    coords_list = [coords] if isinstance(coords, dict) else list(coords)
    cands: list[tuple[list[int], list[str]]] = []
    for c in coords_list:
        cands.extend(_candidates(mol, rings, c, fused))
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
    hetero_by_z = defaultdict(list)  # (b) 按 F>Tl 顺序逐元素收窄该元素原子位次
    for a in all_heteros:
        hetero_by_z[mol.GetAtomWithIdx(a).GetAtomicNum()].append(a)
    for z in _P145_SENIOR:
        if z in hetero_by_z:
            cands = _keep(cands, sorted(hetero_by_z[z]))
    if fused_carbons:
        cands = _keep(cands, fused_carbons)  # (c) 低位次给稠合碳
    if fused_heteros:
        cands = _keep(cands, fused_heteros)  # (d) 低位次给稠合杂原子
    for layer in (sub_layers or ()):
        if len(cands) <= 1:
            break
        if layer:
            cands = _keep(cands, sorted(layer))  # 镜像平局: 逐层按位次集合最小化收窄
    if len(cands) > 1 and alpha_subs:  # P-14.5: 位次集合仍相同时，字母序最前的取代基得最低位次
        def _alpha_key(c):
            """字母序键：[(取代基字母序键, 其位次键)] 排序元组。"""
            return tuple(sorted((k, _label_key(c[1][c[0].index(a)]))
                                for k, a in alpha_subs if a in c[0]))
        best = min(_alpha_key(c) for c in cands)
        cands = [c for c in cands if _alpha_key(c) == best]
    return cands[0]
