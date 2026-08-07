"""Retained chromen-2-one (coumarin) parent (IUPAC P-25 / P-22.2.1 / P-65.6.3).

Fused aromatic 6+6: benzene + α-pyrone lactone (9C + ring O + exocyclic =O).
Numbering fixed: O=1, carbonyl C=2, … 4a, 5–8, 8a. Simple ring prefixes only.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n
from namepredict.layer2.scaffold.naphthalene import _bridge_adjacent, _bridge_pair
from namepredict.layer2.scaffold.quinoline import (
    _bridge_nb_of,
    _path_of_len,
)
from namepredict.layer2.scaffold.ring_parent import (
    _dbl_o_idx,
    _hetero_or_ring_halo,
    _ring_halo_n,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_primary_amines,
    _ring_side_starts,
)
from namepredict.tools.side_alkyl import _side_atoms


def _six_rings(info: dict) -> list[list[int]]:
    rings = info.get("rings") or []
    return [list(r["atom_ids"]) for r in rings if len(r["atom_ids"]) == 6]


def _fused_pair_from(six: list[list[int]], mol: Mol) -> tuple[list[int], list[int], tuple[int, int]] | None:
    for i, r1 in enumerate(six):
        for r2 in six[i + 1 :]:
            bridge = _bridge_pair(r1, r2)
            if bridge is not None and _bridge_adjacent(mol, *bridge):
                return r1, r2, bridge
    return None


def _core_atoms(r1: list[int], r2: list[int]) -> set[int] | None:
    atoms = set(r1) | set(r2)
    return atoms if len(atoms) == 10 else None


def _all_aromatic(mol: Mol, atoms: set[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)


def _one_o_nine_c(mol: Mol, atoms: set[int]) -> int | None:
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms]
    if zs.count(8) != 1 or zs.count(6) != 9:
        return None
    return next(i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)


def _lactone_on_core(info: dict, atoms: set[int]) -> tuple[int, int, int] | None:
    """Return (carbonyl_c, ring_o, a8a) when mono cyclic ester sits on core."""
    esters = info.get("esters") or []
    if len(esters) != 1:
        return None
    e = esters[0]
    c, o, a8a = e["c_idx"], e["o_idx"], e.get("alkoxy_c_idx")
    if c in atoms and o in atoms and a8a in atoms:
        return c, o, a8a
    return None


def _o_next_to_co(mol: Mol, o: int, co: int) -> bool:
    return mol.GetBondBetweenAtoms(o, co) is not None


def _core_from_fused(
    info: dict, r1: list[int], r2: list[int], bridge: tuple[int, int],
) -> tuple[list[int], list[int], int, int, int, int] | None:
    mol, atoms = info["mol"], _core_atoms(r1, r2)
    if atoms is None or not _all_aromatic(mol, atoms):
        return None
    o_idx = _one_o_nine_c(mol, atoms)
    got = _lactone_on_core(info, atoms) if o_idx is not None else None
    if got is None or got[1] != o_idx or not _o_next_to_co(mol, o_idx, got[0]):
        return None
    return r1, r2, o_idx, got[0], bridge[0], bridge[1]


def _chrom_core(info: dict) -> tuple[list[int], list[int], int, int, int, int] | None:
    """Return (r1, r2, o_idx, co_c, ba, bb) or None."""
    fused = _fused_pair_from(_six_rings(info), info["mol"])
    return None if fused is None else _core_from_fused(info, *fused)
