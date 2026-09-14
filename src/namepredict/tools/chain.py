"""碳拓扑原语（tools 层；供 L2 链母体与 L3 侧链共用）。
自 layer2 抽出，供 L3 复用最长开链游走，避免跨层导入。"""
from __future__ import annotations

from rdkit.Chem import Mol


def _carbon_neighbors(mol: Mol, idx: int, banned: set[int] = frozenset()) -> list[int]:
    """链游走用的开链碳邻居，排除 banned 禁走碳。"""
    atom = mol.GetAtomWithIdx(idx)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6 and not n.GetIsAromatic() and not n.IsInRing()
        and n.GetIdx() not in banned
    ]


def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int],
              banned: set[int] = frozenset()) -> list[int]:
    """深度优先寻找从 node 出发的最长开链路径。"""
    best = path
    for nb in _carbon_neighbors(mol, node, banned):
        if nb in path or nb in forbid:
            continue
        cand = _dfs_path(mol, nb, path + [nb], forbid, banned)
        if len(cand) > len(best):
            best = cand
    return best


def _longest_from(mol: Mol, start: int, banned: set[int] = frozenset()) -> list[int]:
    """从 start 出发避开 banned 的最长开链路径。"""
    return _dfs_path(mol, start, [start], set(), banned)
