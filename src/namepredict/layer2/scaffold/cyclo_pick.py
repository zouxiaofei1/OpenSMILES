"""Multi-ring cycloalkane parent selection (L2; P-22.1.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import (
    _all_carbons_are_c,
    _ring_bonds_single,
)

def _sat_c_ring_ok(mol: Mol, atom_ids) -> bool:
    return _all_carbons_are_c(mol, atom_ids) and _ring_bonds_single(mol, atom_ids)


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


