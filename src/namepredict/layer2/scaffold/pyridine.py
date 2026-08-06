"""Retained pyridine / pyridinecarboxylic / pyridinamine / pyridinol / carbonitrile.

Unfused aromatic C5N may be selected among multi-ring mols (P-22.2.1) and may
carry depth-1 aryl arms (Ph / OPh / CH2Ph / OCH2Ph).
"""
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
from namepredict.layer2.aryl_sub import (
    _aryl_atoms,
    _aryl_exclude,
    _aryl_sub_n,
)
from namepredict.layer2.scaffold.ring_parent import (
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


def _is_arom_c5n(mol: Mol, atoms) -> bool:
    if len(atoms) != 6:
        return False
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms]
    if zs.count(7) != 1 or zs.count(6) != 5:
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)


def _is_unfused_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _pyridine_ring_lists(mol: Mol) -> list[list[int]]:
    return [
        list(r) for r in mol.GetRingInfo().AtomRings()
        if _is_arom_c5n(mol, r) and _is_unfused_ring(mol, set(r))
    ]


def _n_in_ring(mol: Mol, ring: list[int]) -> int | None:
    for i in ring:
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 7:
            return i
    return None


def _is_pyridine_core(info: dict) -> bool:
    return bool(_pyridine_ring_lists(info["mol"]))


def _pyridine_n_idx(info: dict, ring: list[int] | None = None) -> int | None:
    if ring is None:
        rings = _pyridine_ring_lists(info["mol"])
        ring = rings[0] if rings else None
    return None if ring is None else _n_in_ring(info["mol"], ring)


def _pyridine_fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_amine"):
        return True
    return bool(info.get("has_ester") or info.get("has_amide") or info.get("has_nitrile"))


def _py_allowed(info: dict, ring_set: set[int]) -> set[int]:
    alk, _ = _arene_alkoxy(info, ring_set)
    return _ring_nitro_atoms(info, ring_set) | alk | _aryl_atoms(info, ring_set)


def _py_ring_ok(info: dict, ring_set: set[int]) -> bool:
    mol: Mol = info["mol"]
    alk, n_alk = _arene_alkoxy(info, ring_set)
    excl = alk | _aryl_exclude(info, ring_set)
    if not _outside_ok(mol, ring_set, _py_allowed(info, ring_set)):
        return False
    return _benzene_subs_ok(
        mol, ring_set, _ring_nitro_n(info, ring_set), n_alk, excl,
        _aryl_sub_n(info, ring_set),
    )


def _pick_pyridine_ring(info: dict) -> list[int] | None:
    mol: Mol = info["mol"]
    cands = [r for r in _pyridine_ring_lists(mol) if _py_ring_ok(info, set(r))]
    if not cands:
        return None
    return max(cands, key=lambda r: (_ring_halo_n(mol, set(r)), -min(r)))


def _is_simple_pyridine(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_block(info):
        return False
    return _pick_pyridine_ring(info) is not None


def _pyridine_parent(info: dict) -> dict:
    ring = _pick_pyridine_ring(info) or _pyridine_ring_lists(info["mol"])[0]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridine",
        "n_idx": _n_in_ring(info["mol"], ring),
    }


def _try_pyridine_parent(info: dict) -> dict | None:
    return _pyridine_parent(info) if _is_simple_pyridine(info) else None


def _pick_one_py_ring(info: dict) -> list[int] | None:
    rings = _pyridine_ring_lists(info["mol"])
    return rings[0] if len(rings) == 1 else None


def _py_cooh_ok(info: dict, ring: list[int]) -> bool:
    mol, ring_set = info["mol"], set(ring)
    if _carboxyl_ring_c(info, ring_set) is None:
        return False
    fg_c = info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring_set, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))


def _is_simple_pyridinecarboxylic(info: dict) -> bool:
    if not _is_pyridine_core(info) or _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    ring = _pick_one_py_ring(info)
    return ring is not None and _py_cooh_ok(info, ring)


def _pyridinecarboxylic_parent(info: dict) -> dict:
    ring = _pick_one_py_ring(info) or list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinecarboxylic",
        "n_idx": _n_in_ring(info["mol"], ring),
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


def _py_fg_full(info: dict, ring_set: set[int], allowed: set[int], alk: set[int]) -> set[int]:
    return allowed | _ring_nitro_atoms(info, ring_set) | alk | _aryl_atoms(info, ring_set)


def _pyridine_subs_ok(mol: Mol, info: dict, ring_set: set[int], allowed: set[int]) -> bool:
    """Allow FG heteroatoms plus aryl/alkoxy/halo/methyl within caps."""
    alk, n_alk = _arene_alkoxy(info, ring_set)
    excl = alk | _aryl_exclude(info, ring_set)
    full = _py_fg_full(info, ring_set, allowed, alk)
    if not _arene_fg_subs_ok(
        mol, ring_set, full, _ring_nitro_n(info, ring_set), 0, n_alk, excl,
        _aryl_sub_n(info, ring_set),
    ):
        return False
    return _ring_halo_n(mol, ring_set) + len(_ring_side_starts(mol, ring_set, excl)) <= 1


def _is_simple_pyridinamine(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_conflict(info, "has_alcohol"):
        return False
    ring = _pick_one_py_ring(info)
    if ring is None:
        return False
    mol, ring_set = info["mol"], set(ring)
    am = _mono_amine_on_ring(info, ring_set)
    if am is None or am.get("degree") != 1:
        return False
    return _pyridine_subs_ok(mol, info, ring_set, {am["n_idx"]})


def _pyridinamine_parent(info: dict) -> dict:
    ring = _pick_one_py_ring(info) or list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinamine",
        "n_idx": _n_in_ring(info["mol"], ring),
        "amine_c_idx": info["amines"][0]["c_idx"],
    }


def _try_pyridinamine_parent(info: dict) -> dict | None:
    return _pyridinamine_parent(info) if _is_simple_pyridinamine(info) else None


def _is_simple_pyridinol(info: dict) -> bool:
    if not _is_pyridine_core(info) or _pyridine_fg_conflict(info, "has_amine"):
        return False
    ring = _pick_one_py_ring(info)
    if ring is None:
        return False
    mol, ring_set = info["mol"], set(ring)
    oh = _mono_oh_on_ring(info, ring_set)
    if oh is None:
        return False
    return _pyridine_subs_ok(mol, info, ring_set, {oh["o_idx"]})


def _pyridinol_parent(info: dict) -> dict:
    ring = _pick_one_py_ring(info) or list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinol",
        "n_idx": _n_in_ring(info["mol"], ring),
        "oh_c_idx": info["hydroxyls"][0]["c_idx"],
    }


def _try_pyridinol_parent(info: dict) -> dict | None:
    return _pyridinol_parent(info) if _is_simple_pyridinol(info) else None


def _try_pyridin_fg_parent(info: dict) -> dict | None:
    """Prefer pyridinamine, else pyridinol (FG-priority parents)."""
    return _try_pyridinamine_parent(info) or _try_pyridinol_parent(info)


def _py_cn_conflict(info: dict) -> bool:
    allow = frozenset({"has_nitrile"})
    return _arene_fg_conflict(info, "has_acid", "has_aldehyde", "has_ketone", allow=allow)


def _py_cn_parts(info: dict, ring: list[int]) -> tuple | None:
    mol, ring_set = info["mol"], set(ring)
    if _fg_ring_c(info, ring_set, "nitriles") is None:
        return None
    fg_c = info["nitriles"][0]["c_idx"]
    n_idx = _nitrile_n_idx(mol, fg_c)
    return (mol, ring_set, fg_c, n_idx, ring) if n_idx is not None else None


def _pyridine_cn_ctx(info: dict) -> tuple | None:
    if not _is_pyridine_core(info) or _py_cn_conflict(info):
        return None
    ring = _pick_one_py_ring(info)
    return None if ring is None else _py_cn_parts(info, ring)


def _is_simple_pyridinecarbonitrile(info: dict) -> bool:
    got = _pyridine_cn_ctx(info)
    if got is None:
        return False
    mol, ring, fg_c, n_idx, _ = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {n_idx})


def _pyridinecarbonitrile_parent(info: dict) -> dict:
    got = _pyridine_cn_ctx(info)
    ring = list(got[4]) if got else list(info["rings"][0]["atom_ids"])
    c = info["nitriles"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyridinecarbonitrile",
        "n_idx": _n_in_ring(info["mol"], ring), "nitrile_c_idx": c,
        "ring_attach_idx": _fg_ring_c(info, set(ring), "nitriles"),
    }


def _try_pyridinecarbonitrile_parent(info: dict) -> dict | None:
    return (
        _pyridinecarbonitrile_parent(info)
        if _is_simple_pyridinecarbonitrile(info) else None
    )
