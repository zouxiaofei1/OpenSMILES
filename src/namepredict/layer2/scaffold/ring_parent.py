from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import Br, C, Cl, F, H, I, O

from namepredict.tools.aryl_sub import (
    _arom_c6_ring_lists, _is_unfused_benzene_ring,
)
from namepredict.tools.side_alkyl import (
    _disjoint_cover, _is_cf3_carbon, _is_cf3_fluoro, _is_omega_halo_c,
    _is_side_halo, _side_sets,
)

def _all_carbons_are_c(mol: Mol, atom_ids: tuple) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == C for i in atom_ids)

def _bond_between(mol: Mol, a: int, b: int):
    return mol.GetBondBetweenAtoms(a, b)

def _ring_bonds_single(mol: Mol, atom_ids: tuple) -> bool:
    ids = list(atom_ids)
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is None or bond.GetBondType().name != "SINGLE":
            return False
    return True

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

def _pure_c_bonds(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        z = other.GetAtomicNum()
        if z not in (1, 6) or (z != 1 and bond.GetBondType().name != "SINGLE"):
            return False
    return True

def _pure_alkyl_outside(mol: Mol, outside: list[int]) -> bool:
    return all(
        _is_cf3_carbon(mol, i) or _is_omega_halo_c(mol, i) or _pure_c_bonds(mol, i)
        for i in outside
    )

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

def _ring_nitro_n(info: dict, ring_set: set[int]) -> int:
    return sum(1 for n in (info.get("nitros") or []) if n["c_idx"] in ring_set)

def _ring_nitro_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for n in info.get("nitros") or []:
        if n["c_idx"] in ring_set:
            out.add(n["n_idx"])
            out.update(n.get("o_idxs") or [])
    return out

def _outside_hetero_ok(atom, ring_set: set[int], allow: set[int]) -> bool:
    if atom.GetAtomicNum() in (1, 6) or atom.GetIdx() in ring_set:
        return True
    if atom.GetIdx() in allow or _is_ring_halo(atom, ring_set):
        return True
    return _is_cf3_fluoro(atom) or _is_side_halo(atom)

def _no_hetero_outside(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    allow = allowed or set()
    return all(_outside_hetero_ok(a, ring_set, allow) for a in mol.GetAtoms())

def _outside_ok(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    if not _no_hetero_outside(mol, ring_set, allowed):
        return False
    return _pure_alkyl_outside(mol, _outside_carbons(mol, ring_set, allowed or set()))

def _ring_double_count(mol: Mol, atom_ids: tuple) -> int:
    ids = list(atom_ids)
    n = 0
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is not None and bond.GetBondType().name == "DOUBLE":
            n += 1
    return n

def _is_carbocycle_ring(info: dict) -> tuple | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    mol: Mol = info["mol"]
    atom_ids = rings[0]["atom_ids"]
    if not _all_carbons_are_c(mol, atom_ids):
        return None
    return atom_ids

def _is_cyclopolyene_core(info: dict) -> bool:
    atom_ids = _is_carbocycle_ring(info)
    if atom_ids is None or info.get("triple_bonds"):
        return False
    return _ring_double_count(info["mol"], atom_ids) >= 2

def _mono_oh_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 1:
        return None
    oh = hydroxyls[0]
    if oh["c_idx"] not in ring_set:
        return None
    return oh

def _mono_amine_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    amines = info.get("amines") or []
    if len(amines) != 1:
        return None
    am = amines[0]
    c_idx = am.get("c_idx")
    if c_idx is None or c_idx not in ring_set:
        return None
    return am

def _dbl_o_idx(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != "DOUBLE":
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == O:
            return other.GetIdx()
    return None

def _is_benzene_core(info: dict) -> bool:
    """True if mol has at least one unfused aromatic C6 carbocycle."""
    mol: Mol = info["mol"]
    for ring in _arom_c6_ring_lists(mol):
        if _is_unfused_benzene_ring(mol, set(ring)):
            return True
    return False

def _ring_halo_n(mol: Mol, ring_set: set[int]) -> int:
    return sum(1 for a in mol.GetAtoms() if _is_ring_halo(a, ring_set))

def _hetero_or_ring_halo(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    return all(_outside_hetero_ok(a, ring_set, allowed) for a in mol.GetAtoms())

def _ring_primary_amines(info: dict, ring_set: set[int]) -> list:
    return [
        a for a in (info.get("amines") or [])
        if a.get("degree") == 1 and a.get("c_idx") in ring_set
    ]

def _is_methyl_on_ring(mol: Mol, s: int, ring_set: set[int]) -> bool:
    if _is_cf3_carbon(mol, s):
        return True
    atom = mol.GetAtomWithIdx(s)
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or (z == 6 and n.GetIdx() in ring_set):
            continue
        return False
    return True

