"""Naphthalene-ring topology primitives (pure rdkit; tools layer).

Extracted from layer2/scaffold/naphthalene.py so L3 heteroaryl side-chain
detection (heteroaryl_sub) can reuse the ring-bridge helpers without importing
a pipeline layer. Only depends on rdkit.
"""
from __future__ import annotations

from rdkit.Chem import Mol


def _bridge_pair(r1: list[int], r2: list[int]) -> tuple[int, int] | None:
    inter = set(r1) & set(r2)
    if len(inter) != 2:
        return None
    a, b = tuple(inter)
    return a, b


def _all_aromatic_c(mol: Mol, atoms: set[int]) -> bool:
    for i in atoms:
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() != 6 or not a.GetIsAromatic():
            return False
    return True


def _bridge_adjacent(mol: Mol, ba: int, bb: int) -> bool:
    return mol.GetBondBetweenAtoms(ba, bb) is not None


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


def _exterior(ring: list[int], ba: int, bb: int) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _path_along(base, ba, bb)
        if p is not None and len(p) == 4:
            return p
    return None


def _exteriors(r1: list[int], r2: list[int], ba: int, bb: int) -> tuple[list[int], list[int]] | None:
    e1, e2 = _exterior(r1, ba, bb), _exterior(r2, ba, bb)
    if e1 is None or e2 is None:
        return None
    return e1, e2


def _naph_chain_from(ba: int, bb: int, e1: list[int], e2: list[int]) -> list[int]:
    """Standard order: 1..4 = e1, 4a=bb, 5..8 = rev(e2), 8a=ba."""
    return e1 + [bb] + list(reversed(e2)) + [ba]


def _chains_for_bridge(
    r1: list[int], r2: list[int], ba: int, bb: int,
) -> list[list[int]]:
    ext = _exteriors(r1, r2, ba, bb)
    if ext is None:
        return []
    e1, e2 = ext
    return [
        _naph_chain_from(ba, bb, e1, e2),
        _naph_chain_from(ba, bb, e2, e1),
    ]
