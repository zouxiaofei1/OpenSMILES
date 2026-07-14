"""Retained 1H-indole parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]pyrrole. NH = 1; mono-methyl / mono-halo only.
Mono ring-C COOH → indolecarboxylic (P-65.1.1).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _arene_subs_ok,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
)
from namepredict.layer2.fused56 import (
    _all_aromatic,
    _chain_atoms,
    _fused56,
    _six_all_c,
)
from namepredict.layer2.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_one_nh(mol: Mol, five: list[int]) -> int | None:
    """Exactly one N (NH, TotalNumHs≥1) and four C on the five-ring."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(7) != 1 or zs.count(6) != 4:
        return None
    n = next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    return n if mol.GetAtomWithIdx(n).GetTotalNumHs() >= 1 else None


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> int | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_one_nh(mol, five)


def _indole_parts(info: dict) -> tuple[list[int], list[int], int, int, int] | None:
    """Return (five, six, nh, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    nh = _core_ok(info["mol"], five, six)
    return None if nh is None else (five, six, nh, ba, bb)


def _is_indole_core(info: dict) -> bool:
    return _indole_parts(info) is not None


def _indole_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_ether", "has_acyl_chloride",
        "has_anhydride",
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


def _indole_subs_ok(mol: Mol, ring: set[int]) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 1:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _is_simple_indole(info: dict) -> bool:
    if not _is_indole_core(info) or _indole_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _indole_parts(info)
    assert parts is not None
    ring = set(parts[0]) | set(parts[1])
    if not _outside_ok(mol, ring):
        return False
    return _indole_subs_ok(mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], nh: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=N, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, nh, ba, bb)


def _indole_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _indole_parts(info)
    assert parts is not None
    five, six, nh, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, nh, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "nh_idx": nh, "bridge": [ba, bb], **extra,
    }


def _indole_parent(info: dict) -> dict:
    return _indole_parent_dict(info, "indole")


def _try_indole_parent(info: dict) -> dict | None:
    return _indole_parent(info) if _is_simple_indole(info) else None


def _indole_ring(info: dict) -> set[int] | None:
    parts = _indole_parts(info)
    return None if parts is None else set(parts[0]) | set(parts[1])


def _is_simple_indolecarboxylic(info: dict) -> bool:
    ring = _indole_ring(info)
    if ring is None or _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    if _carboxyl_ring_c(info, ring) is None:
        return False
    mol, fg_c = info["mol"], info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))


def _try_indolecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_indolecarboxylic(info):
        return None
    ring = _indole_ring(info) or set()
    return _indole_parent_dict(
        info, "indolecarboxylic",
        cooh_c_idx=info["carboxyls"][0]["c_idx"],
        ring_attach_idx=_carboxyl_ring_c(info, ring),
    )
