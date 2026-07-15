"""L2 simple arylboronic acid parent (Ar–B(OH)2, P-68.1).

Gates: mono B(OH)2, c_attach on unfused benzene, outer only simple
prefixes (halo / methyl / CF3; total ≤3). Carboxylic acids outrank.
"""
from __future__ import annotations

from namepredict.layer2.aryl_sub import (
    _arom_c6_ring_lists,
    _is_unfused_benzene_ring,
)
from namepredict.layer2.ring_parent import (
    _arene_fg_subs_ok,
    _is_benzene_core,
    _ring_halo_n,
    _ring_side_starts,
)


_BORONIC_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
    "has_amine", "has_alcohol", "has_thiol", "has_nitro",
)


def _mono_boronic(info: dict) -> dict | None:
    es = info.get("boronics") or []
    return es[0] if info.get("has_boronic") and len(es) == 1 else None


def _no_higher_fgs(info: dict) -> bool:
    return not any(info.get(k) for k in _BORONIC_BAD)


def _boronic_atoms(e: dict) -> set[int]:
    return {e["b_idx"], *e.get("o_idxs", [])}


def _ring_of_attach(mol, c_attach: int) -> set[int] | None:
    for r in _arom_c6_ring_lists(mol):
        rs = set(r)
        if c_attach in rs and _is_unfused_benzene_ring(mol, rs):
            return rs
    return None


def _n_simple_subs(mol, ring: set[int]) -> int:
    return _ring_halo_n(mol, ring) + len(_ring_side_starts(mol, ring))


def _simple_outer_ok(info: dict, ring: set[int], e: dict) -> bool:
    """Only halo / methyl / CF3 outside B(OH)2; total subs ≤3."""
    mol, allow = info["mol"], _boronic_atoms(e)
    if not _arene_fg_subs_ok(mol, ring, allow):
        return False
    return _n_simple_subs(mol, ring) <= 3


def _boronic_ring_ok(info: dict, e: dict) -> set[int] | None:
    ring = _ring_of_attach(info["mol"], e["c_attach"])
    if ring is None:
        return None
    return ring if _simple_outer_ok(info, ring, e) else None


def _is_simple_arylboronic(info: dict) -> bool:
    e = _mono_boronic(info)
    if e is None or not _no_higher_fgs(info) or not _is_benzene_core(info):
        return False
    return _boronic_ring_ok(info, e) is not None


def _pack_boronic(ring: set[int], e: dict) -> dict:
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "boronic",
        "b_idx": e["b_idx"], "c_attach": e["c_attach"],
        "o_idxs": list(e.get("o_idxs") or []),
    }


def _boronic_parent(info: dict) -> dict | None:
    e = _mono_boronic(info)
    if e is None:
        return None
    ring = _boronic_ring_ok(info, e)
    return None if ring is None else _pack_boronic(ring, e)


def _try_boronic(info: dict) -> dict | None:
    return _boronic_parent(info) if _is_simple_arylboronic(info) else None
