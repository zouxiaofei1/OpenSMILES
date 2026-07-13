"""Monocyclic cycloalkanecarboxylic acid parent (IUPAC P-65.1.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import _carboxyl_ring_c, _cooh_oxygen_idxs
from namepredict.layer2.ring_parent import (
    _hetero_or_ring_halo,
    _is_cycloalkane_core,
    _outside_carbons,
    _ring_halo_n,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _walk_linear


def _linear_covers(mol: Mol, start: int, ring: set[int], outside: set[int]) -> bool:
    atoms = _walk_linear(mol, start, ring)
    return atoms is not None and set(atoms) == outside


def _cyclo_acid_alkyl_ok(mol: Mol, ring: set[int], fg_c: int, starts: list[int]) -> bool:
    outside = set(_outside_carbons(mol, ring, {fg_c}))
    if not starts:
        return not outside
    if len(starts) != 1:
        return False
    return _linear_covers(mol, starts[0], ring, outside)


def _cyclo_acid_subs_ok(mol: Mol, ring: set[int], fg_c: int) -> bool:
    starts = _ring_side_starts(mol, ring, {fg_c})
    if _ring_halo_n(mol, ring) + len(starts) > 1:
        return False
    return _cyclo_acid_alkyl_ok(mol, ring, fg_c, starts)


def _is_simple_cycloalkanecarboxylic(info: dict) -> bool:
    if not _is_cycloalkane_core(info) or not info.get("has_acid"):
        return False
    mol: Mol = info["mol"]
    ring = set(info["rings"][0]["atom_ids"])
    if _carboxyl_ring_c(info, ring) is None:
        return False
    fg_c = info["carboxyls"][0]["c_idx"]
    if not _hetero_or_ring_halo(mol, ring, _cooh_oxygen_idxs(mol, fg_c)):
        return False
    return _cyclo_acid_subs_ok(mol, ring, fg_c)


def _cycloalkanecarboxylic_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": len(ring), "kind": "cycloalkanecarboxylic",
        "cooh_c_idx": cooh_c, "ring_attach_idx": _carboxyl_ring_c(info, set(ring)),
    }


def _try_cycloalkanecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_cycloalkanecarboxylic(info):
        return None
    return _cycloalkanecarboxylic_parent(info)
