"""Retained benzofuran / benzofuranamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]furan. O = 1; mono-methyl / mono-halo;
mono primary amine → benzofuranamine.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.fused56 import _chain_atoms, _fused56, _six_all_c, _all_aromatic
from namepredict.layer2.ring_parent import (
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_one_o(mol: Mol, five: list[int]) -> int | None:
    """Exactly one O and four C on the five-ring."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(8) != 1 or zs.count(6) != 4:
        return None
    return next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> int | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_one_o(mol, five)


def _bf_parts(info: dict) -> tuple[list[int], list[int], int, int, int] | None:
    """Return (five, six, o_idx, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    o = _core_ok(info["mol"], five, six)
    return None if o is None else (five, six, o, ba, bb)


def _is_bf_core(info: dict) -> bool:
    return _bf_parts(info) is not None


def _bf_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    if len(starts) != 1:
        return False
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    return outside == starts


def _bf_subs_ok(mol: Mol, ring: set[int]) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 1:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _is_simple_benzofuran(info: dict) -> bool:
    if not _is_bf_core(info) or _bf_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _bf_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    if not _outside_ok(mol, ring):
        return False
    return _bf_subs_ok(mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], o: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=O, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, o, ba, bb)


def _bf_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _bf_parts(info)
    assert parts is not None
    five, six, o, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, o, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "o_idx": o, "bridge": [ba, bb], **extra,
    }


def _benzofuran_parent(info: dict) -> dict:
    return _bf_parent_dict(info, "benzofuran")


def _try_benzofuran_parent(info: dict) -> dict | None:
    return _benzofuran_parent(info) if _is_simple_benzofuran(info) else None


def _bf_amine_conflict(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_thiol",
        "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _bf_primary_amine(info: dict, ring: set[int]) -> dict | None:
    am = _mono_amine_on_ring(info, ring)
    return am if am is not None and am.get("degree") == 1 else None


def _is_simple_benzofuranamine(info: dict) -> bool:
    if not _is_bf_core(info) or _bf_amine_conflict(info):
        return False
    parts = _bf_parts(info)
    assert parts is not None
    mol, ring = info["mol"], _ring_set(parts)
    am = _bf_primary_amine(info, ring)
    if am is None or not _outside_ok(mol, ring, {am["n_idx"]}):
        return False
    return _bf_subs_ok(mol, ring)


def _benzofuranamine_parent(info: dict) -> dict:
    return _bf_parent_dict(
        info, "benzofuranamine", amine_c_idx=info["amines"][0]["c_idx"],
    )


def _try_benzofuranamine_parent(info: dict) -> dict | None:
    return (
        _benzofuranamine_parent(info)
        if _is_simple_benzofuranamine(info) else None
    )
