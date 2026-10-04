"""L2 母体选择所用的脂肪族碳链行走。"""
from __future__ import annotations

from rdkit.Chem import Mol

from opensmiles.tools.chain import _carbon_neighbors, _element_neighbors, _longest_from


def _side_count(mol: Mol, chain: list[int]) -> int:
    """统计链上非链内重原子邻居的个数（支链度）。"""
    chain_set = set(chain)
    n = 0
    for c in chain:
        atom = mol.GetAtomWithIdx(c)
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in chain_set:
                n += 1
    return n


def _better(mol: Mol, cand: list[int], best: list[int]) -> bool:
    if not best:
        return True
    if len(cand) != len(best):
        return len(cand) > len(best)
    return _side_count(mol, cand) > _side_count(mol, best)


def _seed_atoms(mol: Mol, z: int = 6, banned: set[int] = frozenset()) -> list[int]:
    """最长链种子降集：开链子图为多原子树时仅用开链叶，否则回退全集。"""
    atoms = [a.GetIdx() for a in mol.GetAtoms()
             if a.GetAtomicNum() == z and a.GetIdx() not in banned]  # 全同元素原子索引
    if not atoms:
        return []
    has_ring_root = False
    leaves: list[int] = []
    for c in atoms:
        atom = mol.GetAtomWithIdx(c)
        nbs = _element_neighbors(mol, c, z, banned)
        if atom.IsInRing() or atom.GetIsAromatic():
            if nbs:
                has_ring_root = True
        elif len(nbs) == 1:
            leaves.append(c)
    return atoms if (has_ring_root or not leaves) else leaves


def _seed_carbons(mol: Mol, banned: set[int] = frozenset()) -> list[int]:
    """最长碳链的种子降集。"""
    return _seed_atoms(mol, 6, banned)


def _longest_chain(mol: Mol, banned: set[int] = frozenset(), z: int = 6) -> list[int]:
    """返回分子中最长同元素链（降集种子等价加速）。"""
    best: list[int] = []
    for c in _seed_atoms(mol, z, banned):
        path = _longest_from(mol, c, banned=banned, z=z)
        if _better(mol, path, best):
            best = path
    return best


def _component_leaves(mol: Mol, neighbor: int, forbid: int, banned: set[int] = frozenset(),
                      z: int = 6) -> tuple[dict, list[int], int]:
    """DFS neighbor 开链同元素组件（禁走 forbid），返回最深叶。"""
    parent: dict = {neighbor: forbid}
    order = [neighbor]
    dist = {neighbor: 1}
    for x in order:
        for y in _element_neighbors(mol, x, z, banned):
            if y == forbid or y in parent:
                continue
            parent[y] = x
            dist[y] = dist[x] + 1
            order.append(y)
    maxd, leaves = 0, []
    for x in order:
        if any(y != forbid and parent.get(y) == x for y in _element_neighbors(mol, x, z, banned)):
            continue
        d = dist[x]
        if d > maxd:
            maxd, leaves = d, [x]
        elif d == maxd:
            leaves.append(x)
    return parent, leaves, maxd


def _component_path(parent: dict, leaf: int, root: int) -> list[int]:
    """沿 parent 表重建 root→leaf 的链（含两端）。"""
    seg = []
    x = leaf
    while x != root:
        seg.append(x)
        x = parent[x]
    return [root] + list(reversed(seg))


def _all_chains_through(mol: Mol, c_idx: int, banned: set[int] = frozenset(),
                        z: int = 6) -> list[list[int]]:
    """返回 c_idx 的等长最长开链，平局臂全枚举（P-44.4/P-45.2）。"""
    neighbors = _element_neighbors(mol, c_idx, z, banned)
    comps = {n: _component_leaves(mol, n, c_idx, banned, z) for n in neighbors}
    if not comps:
        return [[c_idx]]
    chains: set[tuple[int, ...]] = set()
    ns = list(comps)
    if len(ns) == 1:
        parent, leaves, _ = comps[ns[0]]  # 锚点为链端点：链从最深叶子走向锚点（平局全枚举），叶子在链首、锚点收尾。
        for leaf in leaves:
            chains.add(tuple(reversed(_component_path(parent, leaf, c_idx))))
    else:
        pairs: list[tuple[int, int]] = []  # 锚点在链内：两臂取深度和最大的组件组合（平局全枚举）
        best = 0
        for i in range(len(ns)):
            for j in range(i + 1, len(ns)):
                s = comps[ns[i]][2] + comps[ns[j]][2]
                if s > best:
                    best, pairs = s, [(i, j)]
                elif s == best:
                    pairs.append((i, j))
        for i, j in pairs:
            p1, leaves1, _ = comps[ns[i]]
            p2, leaves2, _ = comps[ns[j]]
            for l1 in leaves1:
                for l2 in leaves2:
                    left = list(reversed(_component_path(p1, l1, c_idx)))[:-1]  # 叶→锚点前
                    right = _component_path(p2, l2, c_idx)                      # 锚点→右叶
                    chains.add(tuple(left + right))
    return [list(c) for c in chains]


def _chain_through_two(mol: Mol, a: int, b: int, banned: set[int] = frozenset(),
                       z: int = 6) -> list[int]:
    """返回同时穿过 c1、c2 的最长链（开链子图是森林，路径唯一）。"""
    if a == b:
        return [a]
    parent, _, _ = _component_leaves(mol, a, -1, banned, z)
    return _component_path(parent, b, a) if b in parent else [a]
