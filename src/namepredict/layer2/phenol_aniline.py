"""Retained phenol / aniline parents with depth-1 aryl (P-63.1.4 / P-62.2.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_sub import (
    _arom_c6_ring_lists,
    _aryl_atoms,
    _aryl_exclude,
    _aryl_sub_n,
    _is_unfused_benzene_ring,
)
from namepredict.layer2.ring_parent import (
    _arene_alkoxy,
    _arene_fg_subs_ok,
    _is_benzene_core,
    _mono_amine_on_ring,
    _mono_oh_on_ring,
    _phenol_amines_ok,
    _ring_alkoxy_atoms,
    _ring_nitro_atoms,
    _ring_nitro_n,
)


def _phenol_allowed(info: dict, ring_set: set[int], oh: dict, ring_ams: list) -> set[int]:
    am_n = {a["n_idx"] for a in ring_ams}
    base = {oh["o_idx"]} | _ring_nitro_atoms(info, ring_set) | am_n
    return base | _ring_alkoxy_atoms(info, ring_set) | _aryl_atoms(info, ring_set)


def _phenol_subs_ok(mol: Mol, info: dict, ring_set: set[int], oh: dict, ring_ams: list) -> bool:
    alk, n_alk = _arene_alkoxy(info, ring_set)
    excl = alk | _aryl_exclude(info, ring_set)
    allowed = _phenol_allowed(info, ring_set, oh, ring_ams)
    return _arene_fg_subs_ok(
        mol, ring_set, allowed, _ring_nitro_n(info, ring_set), len(ring_ams), n_alk, excl,
        _aryl_sub_n(info, ring_set),
    )


def _phenol_ring_ok(info: dict, ring_set: set[int]) -> bool:
    oh, ring_ams = _mono_oh_on_ring(info, ring_set), _phenol_amines_ok(info, ring_set)
    if oh is None or ring_ams is None:
        return False
    return _phenol_subs_ok(info["mol"], info, ring_set, oh, ring_ams)


def _pick_phenol_ring(info: dict) -> set[int] | None:
    mol: Mol = info["mol"]
    for r in _arom_c6_ring_lists(mol):
        rs = set(r)
        if _is_unfused_benzene_ring(mol, rs) and _phenol_ring_ok(info, rs):
            return rs
    return None


def _is_simple_phenol(info: dict) -> bool:
    return bool(_is_benzene_core(info) and _pick_phenol_ring(info) is not None)


def _phenol_parent(info: dict) -> dict | None:
    ring = _pick_phenol_ring(info)
    if ring is None:
        return None
    oh = _mono_oh_on_ring(info, ring)
    return {"chain": list(ring), "n_carbons": 6, "kind": "phenol", "oh_c_idx": oh["c_idx"]}


def _aniline_allowed(info: dict, ring_set: set[int], am: dict) -> set[int]:
    alk, _ = _arene_alkoxy(info, ring_set)
    return {am["n_idx"]} | _ring_nitro_atoms(info, ring_set) | alk | _aryl_atoms(info, ring_set)


def _aniline_subs_ok(info: dict, ring_set: set[int], am: dict) -> bool:
    alk, n_alk = _arene_alkoxy(info, ring_set)
    excl = alk | _aryl_exclude(info, ring_set)
    return _arene_fg_subs_ok(
        info["mol"], ring_set, _aniline_allowed(info, ring_set, am),
        _ring_nitro_n(info, ring_set), 0, n_alk, excl, _aryl_sub_n(info, ring_set),
    )


def _aniline_ring_ok(info: dict, ring_set: set[int]) -> bool:
    if info.get("hydroxyls") or len(info.get("amines") or []) != 1:
        return False
    am = _mono_amine_on_ring(info, ring_set)
    if am is None or am.get("degree") != 1:
        return False
    return _aniline_subs_ok(info, ring_set, am)


def _pick_aniline_ring(info: dict) -> set[int] | None:
    mol: Mol = info["mol"]
    for r in _arom_c6_ring_lists(mol):
        rs = set(r)
        if _is_unfused_benzene_ring(mol, rs) and _aniline_ring_ok(info, rs):
            return rs
    return None


def _is_simple_aniline(info: dict) -> bool:
    return bool(_is_benzene_core(info) and _pick_aniline_ring(info) is not None)


def _aniline_parent(info: dict) -> dict | None:
    ring = _pick_aniline_ring(info)
    if ring is None:
        return None
    am = _mono_amine_on_ring(info, ring)
    return {"chain": list(ring), "n_carbons": 6, "kind": "aniline", "amine_c_idx": am["c_idx"]}
