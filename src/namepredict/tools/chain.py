from __future__ import annotations

from rdkit.Chem import Mol


def _element_neighbors(mol: Mol, idx: int, z: int = 6, banned: set[int] = frozenset()) -> list[int]:
    """链游走用的开链同元素邻居（z = 原子序数），排除 banned 禁走原子。"""
    atom = mol.GetAtomWithIdx(idx)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == z and not n.GetIsAromatic() and not n.IsInRing()
        and n.GetIdx() not in banned
    ]


def _carbon_neighbors(mol: Mol, idx: int, banned: set[int] = frozenset()) -> list[int]:
    """链游走用的开链碳邻居，排除 banned 禁走碳。"""
    return _element_neighbors(mol, idx, 6, banned)


def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int],
              banned: set[int] = frozenset(), z: int = 6) -> list[int]:
    """深度优先寻找从 node 出发的最长开链路径（限同元素 z）。"""
    best = path
    for nb in _element_neighbors(mol, node, z, banned):
        if nb in path or nb in forbid:
            continue
        cand = _dfs_path(mol, nb, path + [nb], forbid, banned, z)
        if len(cand) > len(best):
            best = cand
    return best


def _longest_from(mol: Mol, start: int, banned: set[int] = frozenset(), z: int = 6) -> list[int]:
    """从 start 出发避开 banned 的最长开链路径（限同元素 z）。"""
    return _dfs_path(mol, start, [start], set(), banned, z)
