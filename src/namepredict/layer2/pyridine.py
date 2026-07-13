"""Retained pyridine / pyridinecarboxylic parents (IUPAC P-22.2.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _arene_subs_ok,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
)
from namepredict.layer2.ring_parent import (
    _benzene_subs_ok,
    _outside_ok,
)


def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])


def _is_pyridine_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 6:
        return False
    mol: Mol = info["mol"]
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids]
    if zs.count(7) != 1 or zs.count(6) != 5:
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)


def _pyridine_n_idx(info: dict) -> int | None:
    if not _is_pyridine_core(info):
        return None
    mol: Mol = info["mol"]
    for i in info["rings"][0]["atom_ids"]:
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 7:
            return i
    return None


def _pyridine_fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_amine"):
        return True
    return bool(info.get("has_ester") or info.get("has_amide") or info.get("has_nitrile"))


def _is_simple_pyridine(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _benzene_subs_ok(mol, ring)


def _pyridine_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridine",
        "n_idx": _pyridine_n_idx(info),
    }


def _try_pyridine_parent(info: dict) -> dict | None:
    return _pyridine_parent(info) if _is_simple_pyridine(info) else None


def _is_simple_pyridinecarboxylic(info: dict) -> bool:
    if not _is_pyridine_core(info) or _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if _carboxyl_ring_c(info, ring_set) is None:
        return False
    fg_c = info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring_set, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))


def _pyridinecarboxylic_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinecarboxylic",
        "n_idx": _pyridine_n_idx(info),
        "cooh_c_idx": cooh_c,
        "ring_attach_idx": _carboxyl_ring_c(info, set(ring)),
    }


def _try_pyridinecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_pyridinecarboxylic(info):
        return None
    return _pyridinecarboxylic_parent(info)
