"""Saturated monocyclic hetero carboxylic acids (IUPAC P-65.1.1 / P-22.2.2).

Parent = sat_hetero core with exactly one ring-carbon COOH. Hetero = 1;
COOH locant from ring attach. Optional ≤1 simple ring prefix: hydroxy /
mono-methyl / mono-halo. Chinese retained suffix 甲酸.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
)
from namepredict.layer2.ring_parent import (
    _hetero_or_ring_halo,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)
from namepredict.layer2.sat_hetero import (
    _extra_amine,
    _hetero_asym,
    _hetero_pair,
    _kind_of,
    _mono_methyl_only,
    _n_unsub,
    _ring_atoms_if_mono,
    _ring_n_set,
)

# base sat_hetero kind → carboxylic parent kind
_CARBOXY_KIND = {
    "pyrrolidine": "pyrrolidinecarboxylic",
    "piperidine": "piperidinecarboxylic",
    "piperazine": "piperazinecarboxylic",
    "morpholine": "morpholinecarboxylic",
    "oxolane": "oxolanecarboxylic",
    "oxane": "oxanecarboxylic",
    "thiolane": "thiolanecarboxylic",
    "aziridine": "aziridinecarboxylic",
}


def _shcooh_fg_ok(info: dict) -> bool:
    """Block non-acid principal FGs; alcohol is optional ring prefix only."""
    return not _arene_fg_conflict(
        info, "has_aldehyde", "has_ketone", "has_nitrile",
    )


def _shcooh_base_kind(info: dict) -> str | None:
    kind = _kind_of(info)
    return kind if kind in _CARBOXY_KIND else None


def _shcooh_ctx(info: dict) -> tuple[Mol, list[int], set[int], int, str] | None:
    """Return (mol, atom_ids, ring, fg_c, base_kind) when core + one ring COOH."""
    kind = _shcooh_base_kind(info) if _shcooh_fg_ok(info) else None
    atom_ids = _ring_atoms_if_mono(info) if kind and info.get("carboxyls") else None
    if atom_ids is None:
        return None
    mol, ring = info["mol"], set(atom_ids)
    if _carboxyl_ring_c(info, ring) is None:
        return None
    return mol, atom_ids, ring, info["carboxyls"][0]["c_idx"], kind


def _shcooh_amine_ok(info: dict, mol: Mol, atom_ids: list[int]) -> bool:
    ring_ns = _ring_n_set(mol, atom_ids)
    return not _extra_amine(info, ring_ns) and _n_unsub(mol, ring_ns)


def _shcooh_ohs(info: dict, ring: set[int]) -> list[dict]:
    return [h for h in (info.get("hydroxyls") or []) if h["c_idx"] in ring]


def _shcooh_prefix_n(info: dict, mol: Mol, ring: set[int], fg_c: int) -> int | None:
    """Count halo + methyl + hydroxy prefixes; None if invalid side chains."""
    starts = _ring_side_starts(mol, ring, {fg_c})
    if not _mono_methyl_only(mol, ring | {fg_c}, starts):
        return None
    ohs = _shcooh_ohs(info, ring)
    if len(info.get("hydroxyls") or []) != len(ohs):
        return None
    return _ring_halo_n(mol, ring) + len(starts) + len(ohs)


def _shcooh_allowed(info: dict, mol: Mol, ring: set[int], fg_c: int) -> set[int]:
    ohs = {h["o_idx"] for h in _shcooh_ohs(info, ring)}
    return _cooh_oxygen_idxs(mol, fg_c) | ohs | {fg_c}


def _shcooh_subs_ok(info: dict, mol: Mol, ring: set[int], fg_c: int) -> bool:
    """Allow ≤1 simple ring prefix; COOH oxygens + ring OH O are allowed."""
    n = _shcooh_prefix_n(info, mol, ring, fg_c)
    if n is None or n > 1:
        return False
    allow = _shcooh_allowed(info, mol, ring, fg_c)
    if not _hetero_or_ring_halo(mol, ring, allow):
        return False
    return _outside_ok(mol, ring, allow)


def _is_simple_sat_hetero_carboxylic(info: dict) -> bool:
    ctx = _shcooh_ctx(info)
    if ctx is None:
        return False
    mol, atom_ids, ring, fg_c, _ = ctx
    if not _shcooh_amine_ok(info, mol, atom_ids):
        return False
    return _shcooh_subs_ok(info, mol, ring, fg_c)


def _sat_hetero_carboxylic_parent(info: dict) -> dict:
    _, atom_ids, ring, fg_c, base = _shcooh_ctx(info)
    hs = _hetero_pair(info)
    out = {
        "chain": atom_ids, "n_carbons": len(atom_ids), "kind": _CARBOXY_KIND[base],
        "base_kind": base, "hetero_idxs": hs, "hetero_asym": _hetero_asym(info, hs),
        "cooh_c_idx": fg_c, "ring_attach_idx": _carboxyl_ring_c(info, ring),
    }
    if len(hs) == 1:
        out["hetero_idx"] = hs[0]
    return out


def _try_sat_hetero_carboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_sat_hetero_carboxylic(info):
        return None
    return _sat_hetero_carboxylic_parent(info)
