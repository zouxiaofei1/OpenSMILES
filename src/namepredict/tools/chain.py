"""Open-chain walk primitives (tools layer; shared by L2 chain parents).

Extracted from layer2/chain_walk.py so L3 alkoxy/side detection can reuse
the longest-open-chain walk without importing a pipeline layer.
"""
from __future__ import annotations

from rdkit.Chem import Mol


def _carbon_neighbors(mol: Mol, idx: int) -> list[int]:
    """Open (non-aromatic, non-ring) carbon neighbors for chain FG walks.

    Ring atoms must not enter open parents (else 1-cyclohexylethanone → octan-2-one).
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
