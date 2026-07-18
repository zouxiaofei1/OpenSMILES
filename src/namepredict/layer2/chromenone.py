"""Retained chromen-2-one (coumarin) parent (IUPAC P-25 / P-22.2.1 / P-65.6.3).

Fused aromatic 6+6: benzene + α-pyrone lactone (9C + ring O + exocyclic =O).
Numbering fixed: O=1, carbonyl C=2, … 4a, 5–8, 8a. Simple ring prefixes only.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import _arene_alkoxy
from namepredict.layer2.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n
from namepredict.layer2.naphthalene import _bridge_adjacent, _bridge_pair
from namepredict.layer2.quinoline import (
    _bridge_nb_of,
    _path_of_len,
)
from namepredict.layer2.ring_parent import (
    _dbl_o_idx,
    _hetero_or_ring_halo,
    _ring_halo_n,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_primary_amines,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _side_atoms


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


def _chrom_chain(
    mol: Mol, r1: list[int], r2: list[int], o: int, co: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC: 1=O, 2=carbonyl, 3, 4, 4a, 5..8, 8a (8a = bridge next to O)."""
    nbs = _bridge_nb_of(mol, o, {ba, bb})
    if len(nbs) != 1:
        return None
    a8a, a4a = nbs[0], (bb if nbs[0] == ba else ba)
    pyr, ben = (r1, r2) if o in r1 else (r2, r1)
    mid, ext = _path_of_len(pyr, o, a4a, 3), _path_of_len(ben, a4a, a8a, 4)
    if not mid or not ext or mid[0] != co:
        return None
    return [o] + mid + [a4a] + ext + [a8a]


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_amide",
        "has_nitrile", "has_thiol", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys) or len(info.get("esters") or []) != 1


def _phenol_ohs(info: dict, ring: set[int]) -> list[dict]:
    return [h for h in (info.get("hydroxyls") or []) if h["c_idx"] in ring]


def _chrom_prefix_atoms(info: dict, ring: set[int]) -> set[int]:
    ams = _ring_primary_amines(info, ring)
    ohs = _phenol_ohs(info, ring)
    alk, _ = _arene_alkoxy(info, ring)
    return (
        {a["n_idx"] for a in ams}
        | {h["o_idx"] for h in ohs}
        | _ring_nitro_atoms(info, ring)
        | alk
        | _aryl_atoms(info, ring)
    )


def _chrom_allowed(info: dict, mol: Mol, ring: set[int], co_c: int) -> set[int]:
    o_dbl = _dbl_o_idx(mol, co_c)
    base = {o_dbl} if o_dbl is not None else set()
    return base | _chrom_prefix_atoms(info, ring)


def _chrom_exclude(info: dict, ring: set[int], co_c: int) -> set[int]:
    alk, _ = _arene_alkoxy(info, ring)
    return {co_c} | alk | _aryl_exclude(info, ring)


def _outside_c(mol: Mol, ring: set[int], excl: set[int]) -> set[int]:
    return {
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring and a.GetIdx() not in excl
    }


def _claim_sides(mol: Mol, ring: set[int], starts: list[int]) -> set[int] | None:
    seen: set[int] = set()
    for s in starts:
        atoms = _side_atoms(mol, s, ring)
        if atoms is None or (seen & set(atoms)):
            return None
        seen |= set(atoms)
    return seen


def _side_starts_ok(mol: Mol, ring: set[int], starts: list[int], excl: set[int]) -> bool:
    """Every C-side start is claimable alkyl/prenyl; outside C covered."""
    claimed = _claim_sides(mol, ring, starts)
    return claimed is not None and claimed == _outside_c(mol, ring, excl)


def _prefix_n(info: dict, mol: Mol, ring: set[int], starts: list[int]) -> int:
    n_alk = _arene_alkoxy(info, ring)[1]
    n_am = len(_ring_primary_amines(info, ring)) + len(_phenol_ohs(info, ring))
    return (
        _ring_halo_n(mol, ring) + len(starts) + _ring_nitro_n(info, ring)
        + n_am + n_alk + _aryl_sub_n(info, ring)
    )


def _chrom_subs_ok(info: dict, mol: Mol, ring: set[int], co_c: int) -> bool:
    allowed = _chrom_allowed(info, mol, ring, co_c)
    if not _hetero_or_ring_halo(mol, ring, allowed):
        return False
    excl = _chrom_exclude(info, ring, co_c)
    starts = _ring_side_starts(mol, ring, excl)
    return _side_starts_ok(mol, ring, starts, excl) and _prefix_n(info, mol, ring, starts) <= 3


def _is_simple_chromenone(info: dict) -> bool:
    parts = _chrom_core(info)
    if parts is None or _fg_block(info):
        return False
    mol, ring, co_c = info["mol"], _ring_set(parts), parts[3]
    return _chrom_subs_ok(info, mol, ring, co_c)


def _chromenone_parent(info: dict) -> dict:
    parts = _chrom_core(info)
    assert parts is not None
    r1, r2, o, co, ba, bb = parts
    chain = _chrom_chain(info["mol"], r1, r2, o, co, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 10, "kind": "chromenone",
        "scaffold_id": "chromenone", "o_idx": o, "ketone_c_idx": co,
        "bridge": [ba, bb],
    }


def _try_chromenone_parent(info: dict) -> dict | None:
    return _chromenone_parent(info) if _is_simple_chromenone(info) else None
