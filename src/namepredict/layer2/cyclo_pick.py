"""Multi-ring cycloalkane parent selection (L2; P-22.1.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _all_carbons_are_c,
    _disjoint_cover,
    _outside_carbons,
    _outside_ok,
    _ring_bonds_single,
    _ring_side_starts,
    _side_sets,
)

def _sat_c_ring_ok(mol: Mol, atom_ids) -> bool:
    return _all_carbons_are_c(mol, atom_ids) and _ring_bonds_single(mol, atom_ids)


def _is_cycloalkane_core(info: dict) -> bool:
    """True if ≥1 unfused sat carbocycle exists (multi-ring allowed)."""
    mol: Mol = info["mol"]
    for r in info.get("rings") or []:
        ids = r["atom_ids"]
        if 3 <= len(ids) <= 10 and _sat_c_ring_ok(mol, ids):
            return True
    return False


def _cyclo_ring_unfused(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _one_cyclo_parent_cand(mol: Mol, ids) -> set[int] | None:
    if not (3 <= len(ids) <= 10 and _sat_c_ring_ok(mol, ids)):
        return None
    rs = set(ids)
    return rs if _cyclo_ring_unfused(mol, rs) else None


def _cyclo_parent_candidates(info: dict) -> list[set[int]]:
    mol: Mol = info["mol"]
    out: list[set[int]] = []
    for r in info.get("rings") or []:
        hit = _one_cyclo_parent_cand(mol, r["atom_ids"])
        if hit is not None:
            out.append(hit)
    return out


def _cyclo_sides_cover(mol: Mol, ring_set: set[int]) -> bool:
    """Outside C must be fully covered by recognized alkyl side probes."""
    starts = _ring_side_starts(mol, ring_set)
    outside = set(_outside_carbons(mol, ring_set))
    if not starts:
        return not outside
    sets = _side_sets(mol, ring_set, starts)
    return sets is not None and _disjoint_cover(sets, outside)


def _cyclo_ring_score(mol: Mol, ring: set[int]) -> tuple:
    """Prefer larger ring, then more side starts (as parent)."""
    return (len(ring), len(_ring_side_starts(mol, ring)))


def _cyclo_ring_viable(mol: Mol, ring: set[int]) -> bool:
    return _outside_ok(mol, ring) and _cyclo_sides_cover(mol, ring)


def _better_cyclo(mol: Mol, ring: set[int], best: set[int] | None, sc_best) -> tuple:
    if not _cyclo_ring_viable(mol, ring):
        return best, sc_best
    sc = _cyclo_ring_score(mol, ring)
    if best is None or sc > sc_best:
        return ring, sc
    return best, sc_best


def _pick_cycloalkane_ring(info: dict) -> set[int] | None:
    mol: Mol = info["mol"]
    best, sc_best = None, None
    for ring in _cyclo_parent_candidates(info):
        best, sc_best = _better_cyclo(mol, ring, best, sc_best)
    return best


def _is_simple_cycloalkane(info: dict) -> bool:
    return _pick_cycloalkane_ring(info) is not None


