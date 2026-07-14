"""Aliphatic carbon-chain walks for L2 parent selection."""
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


def _all_carbons(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


def _side_count(mol: Mol, chain: list[int]) -> int:
    chain_set = set(chain)
    n = 0
    for c in chain:
        atom = mol.GetAtomWithIdx(c)
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in chain_set:
                n += 1
    return n


def _chain_key(mol: Mol, path: list[int]) -> tuple:
    return (len(path), _side_count(mol, path))


def _better(mol: Mol, cand: list[int], best: list[int]) -> bool:
    return (not best) or _chain_key(mol, cand) > _chain_key(mol, best)


def _best_among(mol: Mol, seeds: list[int]) -> list[int]:
    best: list[int] = []
    for c in seeds:
        path = _longest_from(mol, c)
        if _better(mol, path, best):
            best = path
    return best


def _longest_chain(mol: Mol, seeds: list[int] | None = None) -> list[int]:
    return _best_among(mol, seeds or _all_carbons(mol))


def _arms_from(mol: Mol, center: int) -> list[list[int]]:
    forbid = {center}
    return [_longest_from(mol, nb, forbid) for nb in _carbon_neighbors(mol, center)]


def _join_through(center: int, arms: list[list[int]]) -> list[int]:
    arms = sorted(arms, key=len, reverse=True)
    if not arms:
        return [center]
    if len(arms) == 1:
        return list(reversed(arms[0])) + [center]
    return list(reversed(arms[0])) + [center] + arms[1]


def _chain_through(info: dict, c_idx: int) -> list[int]:
    mol: Mol = info["mol"]
    return _join_through(c_idx, _arms_from(mol, c_idx))


def _bfs_expand(mol: Mol, cur: int, prev: dict, q: list) -> None:
    for nb in _carbon_neighbors(mol, cur):
        if nb not in prev:
            prev[nb] = cur
            q.append(nb)


def _bfs_prev(mol: Mol, start: int, goal: int) -> dict | None:
    prev: dict = {start: None}
    q = [start]
    while q:
        cur = q.pop(0)
        if cur == goal:
            return prev
        _bfs_expand(mol, cur, prev, q)
    return None


def _rebuild_path(prev: dict, end: int) -> list[int]:
    path = [end]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return list(reversed(path))


def _path_between(mol: Mol, a: int, b: int) -> list[int]:
    if a == b:
        return [a]
    prev = _bfs_prev(mol, a, b)
    return _rebuild_path(prev, b) if prev else [a]


def _best_arm_away(mol: Mol, from_c: int, forbid: set[int]) -> list[int]:
    best: list[int] = []
    for nb in _carbon_neighbors(mol, from_c):
        if nb in forbid:
            continue
        path = _longest_from(mol, nb, forbid | {from_c})
        if len(path) > len(best):
            best = path
    return best


def _chain_through_two(mol: Mol, c1: int, c2: int) -> list[int]:
    path = _path_between(mol, c1, c2)
    left = _best_arm_away(mol, path[0], set(path[1:]))
    right = _best_arm_away(mol, path[-1], set(path[:-1]))
    return list(reversed(left)) + path + right


def _chain_through_bond(mol: Mol, c1: int, c2: int) -> list[int]:
    left = _best_arm_away(mol, c1, {c2})
    right = _best_arm_away(mol, c2, {c1})
    return list(reversed(left)) + [c1, c2] + right
