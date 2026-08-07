"""Unsubstituted sat-hetero monocycle as side chain (piperidinyl etc.)."""
from __future__ import annotations

from collections import Counter

from rdkit.Chem import Mol

# size + hetero z-counts → (en_stem, zh_stem, hetero_Z_for_locant1)
# hetero Z used as numbering origin when present
_SIDE = {
    (5, frozenset([(7, 1)])): ("pyrrolidin", "吡咯烷", 7),
    (6, frozenset([(7, 1)])): ("piperidin", "哌啶", 7),
    (6, frozenset([(8, 1)])): ("oxan", "噁烷", 8),
    (5, frozenset([(8, 1)])): ("oxolan", "氧杂环戊烷", 8),
    (6, frozenset([(8, 1), (7, 1)])): ("morpholin", "吗啉", 7),
}

def _ring_unfused(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True

def _all_single(mol: Mol, ring: set[int]) -> bool:
    ids = list(ring)
    for a in ids:
        for b in ids:
            if a >= b:
                continue
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond is not None and bond.GetBondType().name != "SINGLE":
                return False
    return True

def _z_counts(mol: Mol, ring: set[int]) -> frozenset[tuple[int, int]]:
    c = Counter(
        mol.GetAtomWithIdx(i).GetAtomicNum() for i in ring
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 6
    )
    return frozenset(c.items())

def _ring_sat_ok(mol: Mol, ring: set[int]) -> bool:
    if not (5 <= len(ring) <= 6) or not _ring_unfused(mol, ring):
        return False
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
        return False
    return _all_single(mol, ring)

def _sat_hetero_rings(mol: Mol) -> list[tuple[set[int], tuple]]:
    out = []
    for r in mol.GetRingInfo().AtomRings():
        ring = set(r)
        if not _ring_sat_ok(mol, ring):
            continue
        sig = (len(ring), _z_counts(mol, ring))
        if sig in _SIDE:
            out.append((ring, sig))
    return out

def _exo_only_parent(mol: Mol, ring: set[int], attach: int, parent: set[int]) -> bool:
    for i in ring:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetAtomicNum() == 1 or n.GetIdx() in ring:
                continue
            if not (i == attach and n.GetIdx() in parent):
                return False
    return True

def _hetero1(mol: Mol, ring: set[int], z: int) -> int | None:
    hits = [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == z]
    return hits[0] if len(hits) == 1 else None

def _ring_nbrs(mol: Mol, i: int, ring: set[int]) -> list[int]:
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetIdx() in ring
    ]

def _walk(mol: Mol, ring: set[int], start: int, nxt: int) -> list[int]:
    path, prev, cur = [start], start, nxt
    while cur != start:
        path.append(cur)
        opts = [j for j in _ring_nbrs(mol, cur, ring) if j != prev]
        if not opts:
            break
        prev, cur = cur, opts[0]
    return path

def _locant(mol: Mol, ring: set[int], hetero: int, attach: int) -> int:
    nbrs = _ring_nbrs(mol, hetero, ring)
    if len(nbrs) != 2:
        return 1
    a, b = _walk(mol, ring, hetero, nbrs[0]), _walk(mol, ring, hetero, nbrs[1])
    la = a.index(attach) + 1 if attach in a else 99
    lb = b.index(attach) + 1 if attach in b else 99
    return min(la, lb)

def _side_names(mol: Mol, ring: set[int], start: int, sig: tuple):
    en0, zh0, z1 = _SIDE[sig]
    h1 = _hetero1(mol, ring, z1)
    if h1 is None:
        return None
    loc = _locant(mol, ring, h1, start)
    return f"{en0}-{loc}-yl", f"{zh0}-{loc}-基"

def _try_side(mol: Mol, start: int, parent: set[int], ring: set[int], sig: tuple):
    if start not in ring or ring & parent:
        return None
    if not _exo_only_parent(mol, ring, start, parent):
        return None
    names = _side_names(mol, ring, start, sig)
    return None if names is None else (sorted(ring), names[0], names[1])
