"""Shared fused aromatic 5+6 ring helpers (IUPAC P-22.2.1 / P-25)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _fused_56_pair(rings: list[list[int]]) -> tuple[list[int], list[int]] | None:
    fives = [r for r in rings if len(r) == 5]
    sixes = [r for r in rings if len(r) == 6]
    for five in fives:
        for six in sixes:
            if len(set(five) & set(six)) == 2:
                return five, six
    return None


def _two_rings(info: dict) -> tuple[list[int], list[int]] | None:
    """Adjacent ring pair; multi-ring prefers fused 5+6."""
    rings = [list(r["atom_ids"]) for r in (info.get("rings") or [])]
    if len(rings) < 2:
        return None
    if len(rings) == 2:
        return rings[0], rings[1]
    return _fused_56_pair(rings) or (rings[0], rings[1])


def _size_pair(
    r1: list[int], r2: list[int],
) -> tuple[list[int], list[int]] | None:
    """Return (five, six) if sizes are 5 and 6."""
    if len(r1) == 5 and len(r2) == 6:
        return r1, r2
    if len(r1) == 6 and len(r2) == 5:
        return r2, r1
    return None


def _bridge_pair(r1: list[int], r2: list[int]) -> tuple[int, int] | None:
    inter = set(r1) & set(r2)
    if len(inter) != 2:
        return None
    a, b = tuple(inter)
    return a, b


def _bridge_adjacent(mol: Mol, ba: int, bb: int) -> bool:
    return mol.GetBondBetweenAtoms(ba, bb) is not None


def _all_aromatic(mol: Mol, atoms: set[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)


def _six_all_c(mol: Mol, six: list[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in six)


def _fused_bridge(
    five: list[int], six: list[int], mol: Mol,
) -> tuple[int, int] | None:
    bridge = _bridge_pair(five, six)
    if bridge is None or not _bridge_adjacent(mol, *bridge):
        return None
    return bridge


def _fused56(info: dict) -> tuple[list[int], list[int], tuple[int, int]] | None:
    """Return (five, six, bridge) for adjacent 5+6 fusion, else None."""
    pair = _two_rings(info)
    if pair is None:
        return None
    sized = _size_pair(*pair)
    if sized is None:
        return None
    bridge = _fused_bridge(*sized, info["mol"])
    return None if bridge is None else (sized[0], sized[1], bridge)


def _path_along(ring: list[int], start: int, end: int) -> list[int] | None:
    if start not in ring or end not in ring:
        return None
    i, n = ring.index(start), len(ring)
    path: list[int] = []
    for k in range(1, n):
        atom = ring[(i + k) % n]
        if atom == end:
            return path
        path.append(atom)
    return None


def _path_of_len(
    ring: list[int], start: int, end: int, n: int,
) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _path_along(base, start, end)
        if p is not None and len(p) == n:
            return p
    return None


def _exterior6(six: list[int], a3a: int, a7a: int) -> list[int] | None:
    """Four exterior carbons on the six-ring from 3a toward 7a."""
    return _path_of_len(six, a3a, a7a, 4)


def _hetero_neighbors_on(mol: Mol, h: int, five: list[int]) -> list[int]:
    fset = set(five)
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(h).GetNeighbors()
        if n.GetIdx() in fset
    ]


def _split_hetero_nb(
    mol: Mol, h: int, five: list[int], bridge: set[int],
) -> tuple[int, int] | None:
    """Return (pos2, a7a) where a7a is bridge neighbor of heteroatom."""
    nbs = _hetero_neighbors_on(mol, h, five)
    if len(nbs) != 2:
        return None
    a, b = nbs
    if a in bridge and b not in bridge:
        return b, a
    if b in bridge and a not in bridge:
        return a, b
    return None


def _chain_atoms(
    mol: Mol, five: list[int], six: list[int], h: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=hetero, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    split = _split_hetero_nb(mol, h, five, {ba, bb})
    if split is None:
        return None
    pos2, a7a = split
    a3a = bb if a7a == ba else ba
    mid = _path_of_len(five, pos2, a3a, 1)
    ext = _exterior6(six, a3a, a7a) if mid else None
    return None if not mid or not ext else [h, pos2, mid[0], a3a] + ext + [a7a]
