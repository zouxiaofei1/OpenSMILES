"""开链游走原语（tools 层；供 L2 链母体共享）。
从 layer2/chain_walk.py 抽出，使 L3 烷氧基/侧链检测可复用最长开链游走，无需导入流水线层。"""
from __future__ import annotations

from rdkit.Chem import Mol


def _carbon_neighbors(mol: Mol, idx: int) -> list[int]:
    """链官能团游走用的开链（非芳香、非环）碳邻居。

    环原子不得进入开链母体（否则 1-cyclohexylethanone → octan-2-one）。
    """
    atom = mol.GetAtomWithIdx(idx)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6 and not n.GetIsAromatic() and not n.IsInRing()
    ]


def _extend_best(mol: Mol, node: int, path: list[int], forbid: set[int], best: list[int]) -> list[int]:
    for nb in _carbon_neighbors(mol, node):
        if nb in path or nb in forbid:
            continue
        cand = _dfs_path(mol, nb, path + [nb], forbid)
        if len(cand) > len(best):
            best = cand
    return best


def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int]) -> list[int]:
    return _extend_best(mol, node, path, forbid, path)


def _longest_from(mol: Mol, start: int, forbidden: set[int] | None = None) -> list[int]:
    return _dfs_path(mol, start, [start], forbidden or set())
