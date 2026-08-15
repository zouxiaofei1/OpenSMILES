"""L2 母体选择所用的脂肪族碳链行走。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.tools.chain import _carbon_neighbors, _longest_from


def _all_carbons(mol: Mol) -> list[int]:
    """返回分子中所有碳原子索引列表。"""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


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


def _chain_key(mol: Mol, path: list[int]) -> tuple:
    """构造链比较键：(长度, 支链度)。"""
    return (len(path), _side_count(mol, path))


def _better(mol: Mol, cand: list[int], best: list[int]) -> bool:
    """候选链是否优于当前最佳链（无最佳则胜出）。"""
    return (not best) or _chain_key(mol, cand) > _chain_key(mol, best)


def _best_among(mol: Mol, seeds: list[int]) -> list[int]:
    """在若干种子碳中选出最长链。"""
    best: list[int] = []
    for c in seeds:
        path = _longest_from(mol, c)
        if _better(mol, path, best):
            best = path
    return best


def _longest_chain(mol: Mol, seeds: list[int] | None = None) -> list[int]:
    """返回分子中最长碳链（可选种子约束）。"""
    return _best_among(mol, seeds or _all_carbons(mol))


def _arms_from(mol: Mol, center: int) -> list[list[int]]:
    """以 center 为中心，返回各碳邻居出发的最长臂（禁走 center）。"""
    forbid = {center}
    return [_longest_from(mol, nb, forbid) for nb in _carbon_neighbors(mol, center)]


def _join_through(center: int, arms: list[list[int]]) -> list[int]:
    """按臂长排序拼接穿过 center 的最长链。"""
    arms = sorted(arms, key=len, reverse=True)
    if not arms:
        return [center]
    if len(arms) == 1:
        return list(reversed(arms[0])) + [center]
    return list(reversed(arms[0])) + [center] + arms[1]


def _chain_through(info: dict, c_idx: int) -> list[int]:
    """返回穿过给定碳原子的最长开链。"""
    mol: Mol = info["mol"]
    return _join_through(c_idx, _arms_from(mol, c_idx))


def _bfs_expand(mol: Mol, cur: int, prev: dict, q: list) -> None:
    """BFS 扩展当前节点的碳邻居并记录前驱。"""
    for nb in _carbon_neighbors(mol, cur):
        if nb not in prev:
            prev[nb] = cur
            q.append(nb)


def _bfs_prev(mol: Mol, start: int, goal: int) -> dict | None:
    """BFS 求 start 到 goal 的最短路径前驱表。"""
    prev: dict = {start: None}
    q = [start]
    while q:
        cur = q.pop(0)
        if cur == goal:
            return prev
        _bfs_expand(mol, cur, prev, q)
    return None


def _rebuild_path(prev: dict, end: int) -> list[int]:
    """根据前驱表重建 start→end 的路径。"""
    path = [end]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return list(reversed(path))


def _path_between(mol: Mol, a: int, b: int) -> list[int]:
    """返回两碳原子间的最短路径（含两端）。"""
    if a == b:
        return [a]
    prev = _bfs_prev(mol, a, b)
    return _rebuild_path(prev, b) if prev else [a]


def _best_arm_away(mol: Mol, from_c: int, forbid: set[int]) -> list[int]:
    """返回从 from_c 出发避开 forbid 集合的最长臂。"""
    best: list[int] = []
    for nb in _carbon_neighbors(mol, from_c):
        if nb in forbid:
            continue
        path = _longest_from(mol, nb, forbid | {from_c})
        if len(path) > len(best):
            best = path
    return best


def _chain_through_two(mol: Mol, c1: int, c2: int) -> list[int]:
    """返回同时穿过 c1、c2 的最长链。"""
    path = _path_between(mol, c1, c2)
    left = _best_arm_away(mol, path[0], set(path[1:]))
    right = _best_arm_away(mol, path[-1], set(path[:-1]))
    return list(reversed(left)) + path + right
