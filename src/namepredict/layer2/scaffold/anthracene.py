"""Retained anthracene parent (IUPAC P-25).

Linear fused three C6 aromatic rings (14 C). Unsubstituted or mono-Me/halo.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _anth_shape_ok(s: dict) -> bool:
    return (
        s.get("n_rings") == 3
        and s.get("n_atoms") == 14
        and s.get("topology") == "fused"
        and not s.get("hetero_atoms")
        and s.get("is_aromatic_mancude")
        and len(s.get("fusion_edges") or []) == 2
    )


def _anthracene_system(info: dict) -> dict | None:
    for s in info.get("ring_systems") or []:
        if _anth_shape_ok(s):
            return s
    return None


def _is_linear(system: dict) -> bool:
    edges = system.get("fusion_edges") or []
    if len(edges) != 2:
        return False
    a, b = {edges[0][0], edges[0][1]}, {edges[1][0], edges[1][1]}
    return len(a & b) == 1


def _fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_ether", "has_acyl_chloride",
        "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _subs_ok(mol: Mol, ring: set[int]) -> bool:
    """Unsubstituted only (mono-Me/halo deferred until fixed numbering)."""
    return _ring_halo_n(mol, ring) == 0 and not _ring_side_starts(mol, ring)


def _is_simple_anthracene(info: dict) -> bool:
    system = _anthracene_system(info)
    if system is None or not _is_linear(system) or _fg_block(info):
        return False
    mol, ring = info["mol"], set(system["atom_ids"])
    return _outside_ok(mol, ring) and _subs_ok(mol, ring)


def _anthracene_parent(info: dict) -> dict:
    system = _anthracene_system(info) or {}
    chain = list(system.get("atom_ids") or [])
    return {
        "chain": chain, "n_carbons": 14, "kind": "anthracene",
        "ring_atoms": chain,
    }


def _try_anthracene_parent(info: dict) -> dict | None:
    return _anthracene_parent(info) if _is_simple_anthracene(info) else None
