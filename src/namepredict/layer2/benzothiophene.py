"""Retained 1-benzothiophene / benzothiophenol parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]thiophene. S = 1; mono-methyl / mono-halo;
mono ring OH → 1-benzothiophen-n-ol / 苯并[b]噻吩-n-醇.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.fused56 import _chain_atoms, _fused56, _six_all_c, _all_aromatic
from namepredict.layer2.ring_parent import (
    _mono_oh_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_one_s(mol: Mol, five: list[int]) -> int | None:
    """Exactly one S and four C on the five-ring."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(16) != 1 or zs.count(6) != 4:
        return None
    return next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 16)


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> int | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_one_s(mol, five)


def _bt_parts(info: dict) -> tuple[list[int], list[int], int, int, int] | None:
    """Return (five, six, s_idx, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    s = _core_ok(info["mol"], five, six)
    return None if s is None else (five, six, s, ba, bb)


def _is_bt_core(info: dict) -> bool:
    return _bt_parts(info) is not None


def _bt_fg_block(info: dict) -> bool:
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


def _bt_subs_ok(mol: Mol, ring: set[int]) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 1:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _is_simple_benzothiophene(info: dict) -> bool:
    if not _is_bt_core(info) or _bt_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _bt_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    if not _outside_ok(mol, ring):
        return False
    return _bt_subs_ok(mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], s: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=S, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, s, ba, bb)


def _bt_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _bt_parts(info)
    assert parts is not None
    five, six, s, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, s, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "s_idx": s, "bridge": [ba, bb], **extra,
    }


def _benzothiophene_parent(info: dict) -> dict:
    return _bt_parent_dict(info, "benzothiophene")


def _try_benzothiophene_parent(info: dict) -> dict | None:
    return _benzothiophene_parent(info) if _is_simple_benzothiophene(info) else None


def _bt_ol_conflict(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_ester",
        "has_amide", "has_nitrile", "has_amine", "has_thiol",
        "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _is_simple_benzothiophenol(info: dict) -> bool:
    if not _is_bt_core(info) or _bt_ol_conflict(info):
        return False
    parts = _bt_parts(info)
    assert parts is not None
    mol, ring = info["mol"], _ring_set(parts)
    oh = _mono_oh_on_ring(info, ring)
    if oh is None or not _outside_ok(mol, ring, {oh["o_idx"]}):
        return False
    return _bt_subs_ok(mol, ring)


def _benzothiophenol_parent(info: dict) -> dict:
    return _bt_parent_dict(
        info, "benzothiophenol", oh_c_idx=info["hydroxyls"][0]["c_idx"],
    )


def _try_benzothiophenol_parent(info: dict) -> dict | None:
    return (
        _benzothiophenol_parent(info)
        if _is_simple_benzothiophenol(info) else None
    )
