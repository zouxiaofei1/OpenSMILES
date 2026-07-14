"""Retained 1,3-benzoxazole / benzoxazolamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[d]oxazole. O=1, N=3 (1,3 on five-ring);
≤2 halo / methyl / CF3; mono primary amine at C2 → benzoxazolamine
(with ≤1 halo or CF3).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.fused56 import (
    _all_aromatic,
    _chain_atoms,
    _fused56,
    _six_all_c,
)
from namepredict.layer2.heteroarene5 import _ring_nn_dist
from namepredict.layer2.ring_parent import (
    _is_methyl_on_ring,
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_on(mol: Mol, five: list[int]) -> tuple[int, int] | None:
    """Exactly one O, one N, three C; O–N ring dist 2 (1,3). Return (o, n)."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(8) != 1 or zs.count(7) != 1 or zs.count(6) != 3:
        return None
    o = next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
    n = next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    return (o, n) if _ring_nn_dist(five, [o, n]) == 2 else None


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> tuple[int, int] | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_on(mol, five)


def _box_parts(info: dict) -> tuple[list[int], list[int], int, int, int, int] | None:
    """Return (five, six, o, n, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    on = _core_ok(info["mol"], five, six)
    return None if on is None else (five, six, on[0], on[1], ba, bb)


def _is_box_core(info: dict) -> bool:
    return _box_parts(info) is not None


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _box_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _starts_ok(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring) for s in starts)


def _box_subs_ok(mol: Mol, ring: set[int], cap: int = 2) -> bool:
    """Allow ≤cap simple ring subs: halo + methyl + CF3."""
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > cap:
        return False
    return True if not starts else _starts_ok(mol, ring, starts)


def _is_simple_benzoxazole(info: dict) -> bool:
    if not _is_box_core(info) or _box_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _box_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    return _outside_ok(mol, ring) and _box_subs_ok(mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], o: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=O, 2, 3=N, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, o, ba, bb)


def _box_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _box_parts(info)
    assert parts is not None
    five, six, o, n, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, o, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "o_idx": o, "n_idx": n, "nh_idx": o, "bridge": [ba, bb], **extra,
    }


def _benzoxazole_parent(info: dict) -> dict:
    return _box_parent_dict(info, "benzoxazole")


def _try_benzoxazole_parent(info: dict) -> dict | None:
    return _benzoxazole_parent(info) if _is_simple_benzoxazole(info) else None


def _box_amine_conflict(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_thiol",
        "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _pos2_from_parts(mol: Mol, parts: tuple) -> int | None:
    five, _six, o, n, ba, bb = parts
    bridge = {ba, bb}
    nbs = [
        x.GetIdx() for x in mol.GetAtomWithIdx(o).GetNeighbors()
        if x.GetIdx() in set(five) and x.GetIdx() not in bridge
    ]
    if len(nbs) != 1:
        return None
    return nbs[0] if nbs[0] != n else None


def _primary_amine_on_c2(info: dict, ring: set[int], c2: int) -> dict | None:
    am = _mono_amine_on_ring(info, ring)
    if am is None or am.get("degree") != 1:
        return None
    return am if am.get("c_idx") == c2 else None


def _box_amine_ctx(info: dict) -> tuple[Mol, set[int], dict] | None:
    """Core + C2 primary amine context, or None."""
    if not _is_box_core(info) or _box_amine_conflict(info):
        return None
    parts = _box_parts(info)
    assert parts is not None
    mol, ring = info["mol"], _ring_set(parts)
    c2 = _pos2_from_parts(mol, parts)
    am = _primary_amine_on_c2(info, ring, c2) if c2 is not None else None
    return (mol, ring, am) if am is not None else None


def _is_simple_benzoxazolamine(info: dict) -> bool:
    ctx = _box_amine_ctx(info)
    if ctx is None:
        return False
    mol, ring, am = ctx
    if not _outside_ok(mol, ring, {am["n_idx"]}):
        return False
    return _box_subs_ok(mol, ring, cap=1)


def _benzoxazolamine_parent(info: dict) -> dict:
    return _box_parent_dict(
        info, "benzoxazolamine", amine_c_idx=info["amines"][0]["c_idx"],
    )


def _try_benzoxazolamine_parent(info: dict) -> dict | None:
    return (
        _benzoxazolamine_parent(info)
        if _is_simple_benzoxazolamine(info) else None
    )
