"""Retained 1H-indazole parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[c]pyrazole. NH=1, N=2 (adjacent);
≤2 methyl/halo.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.fused56 import (
    _all_aromatic,
    _chain_atoms,
    _fused56,
    _six_all_c,
)
from namepredict.layer2.scaffold.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_two_n(mol: Mol, five: list[int]) -> list[int] | None:
    """Exactly two N and three C on the five-ring; return N indices."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(7) != 2 or zs.count(6) != 3:
        return None
    return [i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]


def _nn_adjacent(mol: Mol, ns: list[int]) -> bool:
    return mol.GetBondBetweenAtoms(ns[0], ns[1]) is not None


def _nh_of(mol: Mol, ns: list[int]) -> int | None:
    hs = [i for i in ns if mol.GetAtomWithIdx(i).GetTotalNumHs() >= 1]
    return hs[0] if len(hs) == 1 else None


def _five_indazole_nh(mol: Mol, five: list[int]) -> int | None:
    """1H-indazole: adjacent N–N with exactly one NH; return nh_idx."""
    ns = _five_two_n(mol, five)
    if ns is None or not _nn_adjacent(mol, ns):
        return None
    return _nh_of(mol, ns)


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> int | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_indazole_nh(mol, five)


def _iz_parts(info: dict) -> tuple[list[int], list[int], int, int, int] | None:
    """Return (five, six, nh, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    nh = _core_ok(info["mol"], five, six)
    return None if nh is None else (five, six, nh, ba, bb)


def _is_iz_core(info: dict) -> bool:
    return _iz_parts(info) is not None


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _is_methyl(mol: Mol, s: int, ring: set[int]) -> bool:
    return all(
        n.GetAtomicNum() == 1 or n.GetIdx() in ring
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
    )


def _methyl_starts_ok(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl(mol, s, ring) for s in starts)


def _iz_subs_ok(mol: Mol, ring: set[int], cap: int = 2) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > cap:
        return False
    return True if not starts else _methyl_starts_ok(mol, ring, starts)


def _iz_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _is_simple_indazole(info: dict) -> bool:
    if not _is_iz_core(info) or _iz_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _iz_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    return _outside_ok(mol, ring) and _iz_subs_ok(mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], nh: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=NH, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, nh, ba, bb)


def _iz_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _iz_parts(info)
    assert parts is not None
    five, six, nh, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, nh, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "scaffold_id": kind, "nh_idx": nh, "bridge": [ba, bb], **extra,
    }


def _indazole_parent(info: dict) -> dict:
    return _iz_parent_dict(info, "indazole")


def _try_indazole_parent(info: dict) -> dict | None:
    return _indazole_parent(info) if _is_simple_indazole(info) else None
