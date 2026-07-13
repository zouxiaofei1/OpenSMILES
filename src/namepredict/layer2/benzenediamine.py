"""Benzene-1,n-diamine parent (P-62.2.1); mirror of simple benzenediol."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _benzene_alkyl_ns,
    _hetero_or_ring_halo,
    _is_benzene_core,
    _ring_halo_n,
    _ring_primary_amines,
    _ring_side_starts,
)


def _di_amine_on_ring(info: dict, ring_set: set[int]) -> list[dict] | None:
    ams = _ring_primary_amines(info, ring_set)
    if len(ams) != 2 or len(info.get("amines") or []) != 2:
        return None
    return ams


def _bda_starts_ok(mol: Mol, ring_set: set[int], allowed: set[int]) -> list | None:
    starts = _ring_side_starts(mol, ring_set, allowed)
    if len(_benzene_alkyl_ns(mol, ring_set, starts, allowed)) != len(starts):
        return None
    return starts


def _bda_subs_ok(mol: Mol, ring_set: set[int], ams: list) -> bool:
    """≤2 simple halo/methyl (or CF3) outside the two NH2 groups."""
    allowed = {a["n_idx"] for a in ams}
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    starts = _bda_starts_ok(mol, ring_set, allowed)
    return starts is not None and _ring_halo_n(mol, ring_set) + len(starts) <= 2


def _is_simple_benzenediamine(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("hydroxyls"):
        return False
    mol, ring_set = info["mol"], set(info["rings"][0]["atom_ids"])
    ams = _di_amine_on_ring(info, ring_set)
    return ams is not None and _bda_subs_ok(mol, ring_set, ams)


def _benzenediamine_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    ams = _di_amine_on_ring(info, set(ring)) or []
    return {
        "chain": ring,
        "n_carbons": len(ring),
        "kind": "benzenediamine",
        "amine_c_idxs": [a["c_idx"] for a in ams],
    }
