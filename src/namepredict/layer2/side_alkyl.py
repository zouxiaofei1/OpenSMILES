"""Mono alkyl side-chain topology: linear C1–C4, isopropyl or tert-butyl (L2/L3 shared)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _c_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]


def _is_pure_alkyl_c(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    return all(n.GetAtomicNum() in (1, 6) for n in atom.GetNeighbors())


def _nb_kind(n: int, prev: int | None, chain: set[int], start: int, cur: int) -> str:
    if n == prev:
        return "skip"
    if n in chain:
        return "ok" if cur == start else "bad"
    return "free"


def _free_neighbors(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> list[int] | None:
    free: list[int] = []
    for n in _c_neighbors(mol, cur):
        kind = _nb_kind(n, prev, chain, start, cur)
        if kind == "bad":
            return None
        if kind == "free":
            free.append(n)
    return free


def _next_atom(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> int | None | bool:
    free = _free_neighbors(mol, cur, prev, chain, start)
    if free is None or len(free) > 1:
        return False
    return free[0] if free else None


def _advance(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> tuple[int, int | None] | None:
    if not _is_pure_alkyl_c(mol, cur):
        return None
    nxt = _next_atom(mol, cur, prev, chain, start)
    if nxt is False:
        return None
    return cur, nxt  # type: ignore[return-value]


def _walk_linear(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    path: list[int] = []
    prev: int | None = None
    cur: int | None = start
    while cur is not None and len(path) < 5:
        step = _advance(mol, cur, prev, chain, start)
        if step is None:
            return None
        path.append(step[0])
        prev, cur = step[0], step[1]
    return path if 1 <= len(path) <= 4 else None


def _is_terminal_methyl(mol: Mol, idx: int, parent: int) -> bool:
    if not _is_pure_alkyl_c(mol, idx):
        return False
    return _c_neighbors(mol, idx) == [parent]


def _is_isopropyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    free = [n for n in _c_neighbors(mol, start) if n not in chain]
    if len(free) != 2:
        return None
    if not all(_is_terminal_methyl(mol, m, start) for m in free):
        return None
    return [start, free[0], free[1]]


def _is_tert_butyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    free = [n for n in _c_neighbors(mol, start) if n not in chain]
    if len(free) != 3:
        return None
    if not all(_is_terminal_methyl(mol, m, start) for m in free):
        return None
    return [start, *free]


def _side_atoms(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    return (
        _walk_linear(mol, start, chain)
        or _is_isopropyl(mol, start, chain)
        or _is_tert_butyl(mol, start, chain)
    )


def _side_covers(
    mol: Mol, start: int, chain: set[int], outside: set[int],
) -> bool:
    atoms = _side_atoms(mol, start, chain)
    return atoms is not None and set(atoms) == outside


def _side_sets(mol: Mol, chain: set[int], starts: list[int]) -> list[set[int]] | None:
    sets: list[set[int]] = []
    for s in starts:
        atoms = _side_atoms(mol, s, chain)
        if atoms is None:
            return None
        sets.append(set(atoms))
    return sets


def _disjoint_cover(sets: list[set[int]], outside: set[int]) -> bool:
    seen: set[int] = set()
    for s in sets:
        if seen & s:
            return False
        seen |= s
    return seen == outside
