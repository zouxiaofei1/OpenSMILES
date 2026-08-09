"""Aryl arms: Ph / OPh / CH2Ph / OCH2Ph with recursive nested leaves (P-29.3).

Ph may carry ≤3 leaves (simple + nested Ph/OPh up to max depth).
Locants from attach=1 (lowest set).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C, H

from namepredict.constants import HALO_Z as _HALO

def _is_arom_c6(mol: Mol, atoms) -> bool:
    if len(atoms) != 6:
        return False
    return all(
        mol.GetAtomWithIdx(i).GetIsAromatic()
        and mol.GetAtomWithIdx(i).GetAtomicNum() == C
        for i in atoms
    )

def _c6_rings_at(mol: Mol, c_idx: int) -> list[set[int]]:
    return [
        set(r) for r in mol.GetRingInfo().AtomRings()
        if c_idx in r and _is_arom_c6(mol, r)
    ]

def _nb_outside(mol: Mol, i: int, ring: set[int]):
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != H and n.GetIdx() not in ring
    ]



def _is_unfused_benzene_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _arom_c6_ring_lists(mol: Mol) -> list[list[int]]:
    return [list(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]

def _unfused_c6_at(mol: Mol, c_idx: int) -> set[int] | None:
    """Sole unfused aromatic C6 containing c_idx, else None."""
    hits = [
        set(r) for r in _arom_c6_ring_lists(mol)
        if c_idx in r and _is_unfused_benzene_ring(mol, set(r))
    ]
    return hits[0] if len(hits) == 1 else None

def _exocyclic_fg_ring(mol: Mol, fg_c: int) -> set[int] | None:
    """Unfused C6 attached to exocyclic FG carbon (COOH/CHO/CN/...)."""
    nbs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(fg_c).GetNeighbors()
        if n.GetAtomicNum() == C and n.GetIsAromatic()
    ]
    return _unfused_c6_at(mol, nbs[0]) if len(nbs) == 1 else None
