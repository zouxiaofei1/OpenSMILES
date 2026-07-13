"""Retained 1H-indole parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]pyrrole. NH = 1; mono-methyl / mono-halo only.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _two_rings(info: dict) -> tuple[list[int], list[int]] | None:
    rings = info.get("rings") or []
    if len(rings) != 2:
        return None
    return list(rings[0]["atom_ids"]), list(rings[1]["atom_ids"])


def _size_pair(
    r1: list[int], r2: list[int],
) -> tuple[list[int], list[int]] | None:
    """Return (five, six) if sizes are 5 and 6."""
    if len(r1) == 5 and len(r2) == 6:
        return r1, r2
    if len(r1) == 6 and len(r2) == 5:
        return r2, r1
    return None


def _bridge_pair(r1: list[int], r2: list[int]) -> tuple[int, int] | None:
    inter = set(r1) & set(r2)
    if len(inter) != 2:
        return None
    a, b = tuple(inter)
    return a, b


def _bridge_adjacent(mol: Mol, ba: int, bb: int) -> bool:
    return mol.GetBondBetweenAtoms(ba, bb) is not None


def _all_aromatic(mol: Mol, atoms: set[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)


def _six_all_c(mol: Mol, six: list[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in six)


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


def _fused_bridge(
    five: list[int], six: list[int], mol: Mol,
) -> tuple[int, int] | None:
    bridge = _bridge_pair(five, six)
    if bridge is None or not _bridge_adjacent(mol, *bridge):
        return None
    return bridge


def _fused56(info: dict) -> tuple[list[int], list[int], tuple[int, int]] | None:
    """Return (five, six, bridge) for adjacent 5+6 fusion, else None."""
    pair = _two_rings(info)
    if pair is None:
        return None
    sized = _size_pair(*pair)
    if sized is None:
        return None
    bridge = _fused_bridge(*sized, info["mol"])
    return None if bridge is None else (sized[0], sized[1], bridge)


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


def _path_along(ring: list[int], start: int, end: int) -> list[int] | None:
    if start not in ring or end not in ring:
        return None
    i, n = ring.index(start), len(ring)
    path: list[int] = []
    for k in range(1, n):
        atom = ring[(i + k) % n]
        if atom == end:
            return path
        path.append(atom)
    return None


def _path_of_len(
    ring: list[int], start: int, end: int, n: int,
) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _path_along(base, start, end)
        if p is not None and len(p) == n:
            return p
    return None


def _exterior6(six: list[int], a3a: int, a7a: int) -> list[int] | None:
    """Four exterior carbons on the six-ring from 3a toward 7a."""
    return _path_of_len(six, a3a, a7a, 4)


def _nh_neighbors_on(mol: Mol, nh: int, five: list[int]) -> list[int]:
    fset = set(five)
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(nh).GetNeighbors()
        if n.GetIdx() in fset
    ]


def _split_nh_nb(mol: Mol, nh: int, five: list[int], bridge: set[int]) -> tuple[int, int] | None:
    """Return (pos2, a7a) where a7a is bridge neighbor of NH."""
    nbs = _nh_neighbors_on(mol, nh, five)
    if len(nbs) != 2:
        return None
    a, b = nbs
    if a in bridge and b not in bridge:
        return b, a
    if b in bridge and a not in bridge:
        return a, b
    return None


def _chain_atoms(
    mol: Mol, five: list[int], six: list[int], nh: int, ba: int, bb: int,
) -> list[int] | None:
    split = _split_nh_nb(mol, nh, five, {ba, bb})
    if split is None:
        return None
    pos2, a7a = split
    a3a = bb if a7a == ba else ba
    mid = _path_of_len(five, pos2, a3a, 1)
    ext = _exterior6(six, a3a, a7a) if mid else None
    return None if not mid or not ext else [nh, pos2, mid[0], a3a] + ext + [a7a]


def _build_chain(
    mol: Mol, five: list[int], six: list[int], nh: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=N, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, nh, ba, bb)


def _indole_parent(info: dict) -> dict:
    parts = _indole_parts(info)
    assert parts is not None
    five, six, nh, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, nh, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": "indole",
        "nh_idx": nh, "bridge": [ba, bb],
    }


def _try_indole_parent(info: dict) -> dict | None:
    return _indole_parent(info) if _is_simple_indole(info) else None
