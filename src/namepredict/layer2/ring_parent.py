from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import Br, C, Cl, F, H, I, O

def _outside_carbons(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    return [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == C
        and a.GetIdx() not in ring_set
        and a.GetIdx() not in skip
    ]

def _ring_side_starts(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    starts: list[int] = []
    for r in ring_set:
        for n in mol.GetAtomWithIdx(r).GetNeighbors():
            if n.GetAtomicNum() == C and n.GetIdx() not in ring_set:
                if n.GetIdx() not in skip:
                    starts.append(n.GetIdx())
    return starts

def _is_ring_halo(atom, ring_set: set[int]) -> bool:
    if atom.GetAtomicNum() not in (F, Cl, Br, I):
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetIdx() in ring_set

def _no_hetero_outside(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    allow = allowed or set()
    return all(_outside_hetero_ok(a, ring_set, allow) for a in mol.GetAtoms())

def _outside_ok(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    if not _no_hetero_outside(mol, ring_set, allowed):
        return False
    return _pure_alkyl_outside(mol, _outside_carbons(mol, ring_set, allowed or set()))

def _o_idx(mol: Mol, c_idx: int, bond_name: str) -> int | None:
    """Index of O neighbour of `c_idx` across a bond of the given type."""
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != bond_name:
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == O:
            return other.GetIdx()
    return None


def _dbl_o_idx(mol: Mol, c_idx: int) -> int | None:
    return _o_idx(mol, c_idx, "DOUBLE")

def _ring_halo_n(mol: Mol, ring_set: set[int]) -> int:
    return sum(1 for a in mol.GetAtoms() if _is_ring_halo(a, ring_set))


