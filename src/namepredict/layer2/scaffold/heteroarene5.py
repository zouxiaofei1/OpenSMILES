"""Retained monocyclic heteroarene parents (IUPAC P-22.2.1 / P-22.1).

Five-membered mono-hetero: furan / thiophene / 1H-pyrrole (optional monomethyl/monohalo).
Five-membered di-aza: 1H-imidazole / 1H-pyrazole (optional monomethyl/monohalo).
Six-membered diazines: pyrimidine / pyrazine / pyridazine (simple ring subs / pyrimidinamine).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import (
    _arene_alkoxy,
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _linear_or_omega_halo_sides_ok
from namepredict.layer2.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n

# Z of ring hetero → parent kind (5-membered mono)
_HETERO5_KIND = {8: "furan", 16: "thiophene", 7: "pyrrole"}
# min ring distance between two N → diazine kind
_DIAZINE_KIND = {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"}

def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])

def _hetero5_zs(mol: Mol, atom_ids: list[int]) -> list[int]:
    return [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids]

def _hetero5_fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    return bool(info.get("has_nitrile") or info.get("has_thiol") or info.get("has_nitro"))

def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:
    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])
    d = abs(ia - ib)
    return min(d, len(atom_ids) - d)

def _is_diazole_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 5:
        return False
    mol: Mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
        return False
    zs = _hetero5_zs(mol, atom_ids)
    return zs.count(7) == 2 and zs.count(6) == 3

def _ring_n_idxs(mol: Mol, ring: list[int]) -> list[int]:
    return [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]

def _nh_among(mol: Mol, ns: list[int]) -> int | None:
    hs = [i for i in ns if mol.GetAtomWithIdx(i).GetTotalNumHs() >= 1]
    return hs[0] if len(hs) == 1 else None

def _diazole_n_pair(info: dict, nn_dist: int) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-diazole with given N–N ring distance."""
    if not _is_diazole_core(info):
        return None
    mol, ring = info["mol"], list(info["rings"][0]["atom_ids"])
    ns = _ring_n_idxs(mol, ring)
    if len(ns) != 2 or _ring_nn_dist(ring, ns) != nn_dist:
        return None
    nh = _nh_among(mol, ns)
    return None if nh is None else (nh, ns[0] if ns[1] == nh else ns[1])

def _pyrazole_n_pair(info: dict) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-pyrazole (1,2-diazole); else None."""
    return _diazole_n_pair(info, 1)
