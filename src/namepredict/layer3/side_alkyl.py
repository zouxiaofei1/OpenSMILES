"""Mono alkyl side-chain topology: linear C1–C4 + retained branched (L2/L3)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C, F, H, HALO_Z as _HALO_Z



def _is_side_halo(atom) -> bool:
    """Terminal halogen on a non-ring carbon (ω-haloalkyl)."""
    if atom.GetAtomicNum() not in _HALO_Z:
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]
    if len(heavies) != 1 or heavies[0].GetAtomicNum() != C:
        return False
    return not heavies[0].IsInRing()

def _first_side(mol: Mol, start: int, chain: set[int], probes) -> list[int] | None:
    for probe in probes:
        got = probe(mol, start, chain)
        if got is not None:
            return got
    return None

def _side_atoms(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    return _first_side(mol, start, chain, ())

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

