"""Retained naphthalene parent (IUPAC P-22.1.1 / P-25).

Two six-membered aromatic carbocycles sharing exactly two adjacent atoms.
Unsubstituted, mono-methyl, or mono-halo only.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _two_six_rings(info: dict) -> tuple[list[int], list[int]] | None:
    rings = info.get("rings") or []
    if len(rings) != 2:
        return None
    r1, r2 = list(rings[0]["atom_ids"]), list(rings[1]["atom_ids"])
    return (r1, r2) if len(r1) == 6 and len(r2) == 6 else None


def _bridge_pair(r1: list[int], r2: list[int]) -> tuple[int, int] | None:
    inter = set(r1) & set(r2)
    if len(inter) != 2:
        return None
    a, b = tuple(inter)
    return a, b


def _all_aromatic_c(mol: Mol, atoms: set[int]) -> bool:
    for i in atoms:
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() != 6 or not a.GetIsAromatic():
            return False
    return True


def _bridge_adjacent(mol: Mol, ba: int, bb: int) -> bool:
    return mol.GetBondBetweenAtoms(ba, bb) is not None


def _core_atoms_ok(mol: Mol, r1: list[int], r2: list[int]) -> bool:
    atoms = set(r1) | set(r2)
    return len(atoms) == 10 and _all_aromatic_c(mol, atoms)


def _is_naphthalene_core(info: dict) -> bool:
    pair = _two_six_rings(info)
    if pair is None:
        return False
    bridge = _bridge_pair(*pair)
    if bridge is None:
        return False
    mol: Mol = info["mol"]
    return _core_atoms_ok(mol, *pair) and _bridge_adjacent(mol, *bridge)


def _naph_fg_block(info: dict) -> bool:
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


def _naph_subs_ok(mol: Mol, ring: set[int]) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 1:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _is_simple_naphthalene(info: dict) -> bool:
    if not _is_naphthalene_core(info) or _naph_fg_block(info):
        return False
    mol: Mol = info["mol"]
    ring = set(info["rings"][0]["atom_ids"]) | set(info["rings"][1]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _naph_subs_ok(mol, ring)


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


def _exterior(ring: list[int], ba: int, bb: int) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _path_along(base, ba, bb)
        if p is not None and len(p) == 4:
            return p
    return None


def _exteriors(r1: list[int], r2: list[int], ba: int, bb: int) -> tuple[list[int], list[int]] | None:
    e1, e2 = _exterior(r1, ba, bb), _exterior(r2, ba, bb)
    if e1 is None or e2 is None:
        return None
    return e1, e2


def _naph_chain_from(ba: int, bb: int, e1: list[int], e2: list[int]) -> list[int]:
    """Standard order: 1..4 = e1, 4a=bb, 5..8 = rev(e2), 8a=ba."""
    return e1 + [bb] + list(reversed(e2)) + [ba]


def _chains_for_bridge(
    r1: list[int], r2: list[int], ba: int, bb: int,
) -> list[list[int]]:
    ext = _exteriors(r1, r2, ba, bb)
    if ext is None:
        return []
    e1, e2 = ext
    return [
        _naph_chain_from(ba, bb, e1, e2),
        _naph_chain_from(ba, bb, e2, e1),
    ]


def _naph_chains(info: dict) -> list[list[int]]:
    pair = _two_six_rings(info)
    if pair is None:
        return []
    bridge = _bridge_pair(*pair)
    if bridge is None:
        return []
    ba, bb = bridge
    out = _chains_for_bridge(*pair, ba, bb)
    out.extend(_chains_for_bridge(*pair, bb, ba))
    return out


def _naphthalene_parent(info: dict) -> dict:
    chains = _naph_chains(info)
    chain = chains[0] if chains else []
    return {
        "chain": chain, "n_carbons": 10, "kind": "naphthalene",
        "bridge": list(_bridge_pair(*_two_six_rings(info)) or ()),
        "naph_chains": chains,
    }


def _try_naphthalene_parent(info: dict) -> dict | None:
    return _naphthalene_parent(info) if _is_simple_naphthalene(info) else None
