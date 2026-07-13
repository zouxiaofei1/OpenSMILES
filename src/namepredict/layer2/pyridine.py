"""Retained pyridine / pyridinecarboxylic / pyridinamine / pyridinol / carbonitrile."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _arene_subs_ok,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
    _fg_ring_c,
    _nitrile_n_idx,
)
from namepredict.layer2.ring_parent import (
    _arene_alkoxy,
    _arene_fg_subs_ok,
    _benzene_subs_ok,
    _mono_amine_on_ring,
    _mono_oh_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_side_starts,
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
    alk, n_alk = _arene_alkoxy(info, ring)
    allowed = _ring_nitro_atoms(info, ring) | alk
    if not _outside_ok(mol, ring, allowed):
        return False
    return _benzene_subs_ok(mol, ring, _ring_nitro_n(info, ring), n_alk, alk)


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


def _pyridine_fg_conflict(info: dict, *extra: str) -> bool:
    bad = (
        "has_acid", "has_ester", "has_amide", "has_nitrile",
        "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
        *extra,
    )
    return any(info.get(k) for k in bad)


def _pyridine_subs_ok(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    """Allow FG heteroatoms plus at most one ring mono-halo or mono-methyl."""
    if not _arene_fg_subs_ok(mol, ring_set, allowed, 0, 0, 0, None):
        return False
    starts = _ring_side_starts(mol, ring_set)
    return _ring_halo_n(mol, ring_set) + len(starts) <= 1


def _is_simple_pyridinamine(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_conflict(info, "has_alcohol"):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    am = _mono_amine_on_ring(info, ring)
    if am is None or am.get("degree") != 1:
        return False
    return _pyridine_subs_ok(mol, ring, {am["n_idx"]})


def _pyridinamine_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinamine",
        "n_idx": _pyridine_n_idx(info),
        "amine_c_idx": info["amines"][0]["c_idx"],
    }


def _try_pyridinamine_parent(info: dict) -> dict | None:
    return _pyridinamine_parent(info) if _is_simple_pyridinamine(info) else None


def _is_simple_pyridinol(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_conflict(info, "has_amine"):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    oh = _mono_oh_on_ring(info, ring)
    if oh is None:
        return False
    return _pyridine_subs_ok(mol, ring, {oh["o_idx"]})


def _pyridinol_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinol",
        "n_idx": _pyridine_n_idx(info),
        "oh_c_idx": info["hydroxyls"][0]["c_idx"],
    }


def _try_pyridinol_parent(info: dict) -> dict | None:
    return _pyridinol_parent(info) if _is_simple_pyridinol(info) else None


def _try_pyridin_fg_parent(info: dict) -> dict | None:
    """Prefer pyridinamine, else pyridinol (FG-priority parents)."""
    return _try_pyridinamine_parent(info) or _try_pyridinol_parent(info)


def _pyridine_cn_ctx(info: dict) -> tuple | None:
    if not _is_pyridine_core(info):
        return None
    allow = frozenset({"has_nitrile"})
    if _arene_fg_conflict(info, "has_acid", "has_aldehyde", "has_ketone", allow=allow):
        return None
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if _fg_ring_c(info, ring, "nitriles") is None:
        return None
    fg_c, n_idx = info["nitriles"][0]["c_idx"], _nitrile_n_idx(info["mol"], info["nitriles"][0]["c_idx"])
    return (mol, ring, fg_c, n_idx) if n_idx is not None else None


def _is_simple_pyridinecarbonitrile(info: dict) -> bool:
    got = _pyridine_cn_ctx(info)
    if got is None:
        return False
    mol, ring, fg_c, n_idx = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {n_idx})


def _pyridinecarbonitrile_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    c = info["nitriles"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinecarbonitrile",
        "n_idx": _pyridine_n_idx(info), "nitrile_c_idx": c,
        "ring_attach_idx": _fg_ring_c(info, set(ring), "nitriles"),
    }


def _try_pyridinecarbonitrile_parent(info: dict) -> dict | None:
    return (
        _pyridinecarbonitrile_parent(info)
        if _is_simple_pyridinecarbonitrile(info) else None
    )
