"""Pick simple benzene parent ring among unfused aromatic C6 (P-22 / P-29.3)."""
from __future__ import annotations

from namepredict.layer2.aryl_sub import (
    _aryl_atoms, _aryl_sub_n, _arom_c6_ring_lists, _is_unfused_benzene_ring,
    _phenyl_starts_set,
)
from namepredict.layer2.scaffold.ring_parent import (
    _arene_alkoxy, _benzene_subs_ok, _outside_ok, _ring_halo_n,
    _ring_iso_atoms, _ring_iso_n, _ring_nitro_atoms, _ring_nitro_n,
    _ring_side_starts,
)


def _benzene_exclude(info: dict, ring_set: set[int]) -> set[int]:
    mol = info["mol"]
    return (
        _ring_nitro_atoms(info, ring_set)
        | _ring_iso_atoms(info, ring_set)
        | _ring_alkoxy_atoms_local(info, ring_set)
        | _aryl_atoms(info, ring_set)
        | _phenyl_starts_set(mol, ring_set)
    )


def _ring_alkoxy_atoms_local(info: dict, ring_set: set[int]) -> set[int]:
    return _arene_alkoxy(info, ring_set)[0]


def _benzene_side_excl(info: dict, ring_set: set[int], alk: set[int]) -> set[int]:
    mol = info["mol"]
    return (
        alk | _aryl_atoms(info, ring_set) | _phenyl_starts_set(mol, ring_set)
        | _ring_iso_atoms(info, ring_set)
    )


def _benzene_ring_ok(info: dict, ring_set: set[int]) -> bool:
    mol = info["mol"]
    alk, n_alk = _arene_alkoxy(info, ring_set)
    if not _outside_ok(mol, ring_set, _benzene_exclude(info, ring_set)):
        return False
    excl = _benzene_side_excl(info, ring_set, alk)
    return _benzene_subs_ok(
        mol, ring_set, _ring_nitro_n(info, ring_set), n_alk, excl,
        _aryl_sub_n(info, ring_set), _ring_iso_n(info, ring_set),
    )


def _benzene_ring_score(info: dict, ring: set[int]) -> tuple:
    mol = info["mol"]
    excl = _aryl_atoms(info, ring) | _phenyl_starts_set(mol, ring)
    starts = _ring_side_starts(mol, ring, excl)
    return (_ring_halo_n(mol, ring), len(starts), -min(ring))


def _pick_benzene_ring(info: dict) -> list[int] | None:
    mol = info["mol"]
    cands = [
        r for r in _arom_c6_ring_lists(mol)
        if _is_unfused_benzene_ring(mol, set(r)) and _benzene_ring_ok(info, set(r))
    ]
    if not cands:
        return None
    return max(cands, key=lambda r: _benzene_ring_score(info, set(r)))


def _is_simple_benzene(info: dict) -> bool:
    return _pick_benzene_ring(info) is not None
