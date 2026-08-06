"""Retained quinoline / isoquinoline parent (IUPAC P-22.2.1 / P-25).

Two six-membered aromatic rings (9C+1N) sharing two adjacent bridge atoms.
N adjacent to one bridge → quinoline (N=1); N adjacent to none → isoquinoline (N=2).
≤2 halo/methyl; mono ring OH → quinolinol; mono COOH → quinolinecarboxylic.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _arene_subs_ok,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
)
from namepredict.layer2.scaffold.naphthalene import (
    _bridge_adjacent,
    _bridge_pair,
    _two_six_rings,
)
from namepredict.layer2.scaffold.ring_parent import (
    _is_methyl_on_ring,
    _mono_oh_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)

def _core_atoms(r1: list[int], r2: list[int]) -> set[int] | None:
    atoms = set(r1) | set(r2)
    return atoms if len(atoms) == 10 else None

def _all_aromatic(mol: Mol, atoms: set[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)

def _one_n_nine_c(mol: Mol, atoms: set[int]) -> int | None:
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms]
    if zs.count(7) != 1 or zs.count(6) != 9:
        return None
    return next(i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)

def _n_bridge_degree(n_idx: int, bridge: set[int], mol: Mol) -> int:
    nbs = {a.GetIdx() for a in mol.GetAtomWithIdx(n_idx).GetNeighbors()}
    return len(nbs & bridge)

def _kind_from_n(n_idx: int, bridge: set[int], mol: Mol) -> str | None:
    if n_idx in bridge:
        return None
    deg = _n_bridge_degree(n_idx, bridge, mol)
    if deg == 1:
        return "quinoline"
    return "isoquinoline" if deg == 0 else None

def _fused_pair(info: dict) -> tuple[list[int], list[int], tuple[int, int]] | None:
    pair = _two_six_rings(info)
    if pair is None:
        return None
    bridge = _bridge_pair(*pair)
    if bridge is None or not _bridge_adjacent(info["mol"], *bridge):
        return None
    return (*pair, bridge)

def _q_from_fused(
    info: dict, r1: list[int], r2: list[int], bridge: tuple[int, int],
) -> tuple[list[int], list[int], int, int, int, str] | None:
    atoms = _core_atoms(r1, r2)
    if atoms is None or not _all_aromatic(info["mol"], atoms):
        return None
    n_idx = _one_n_nine_c(info["mol"], atoms)
    if n_idx is None:
        return None
    kind = _kind_from_n(n_idx, set(bridge), info["mol"])
    return None if kind is None else (r1, r2, n_idx, bridge[0], bridge[1], kind)

def _q_core(info: dict) -> tuple[list[int], list[int], int, int, int, str] | None:
    """Return (r1, r2, n_idx, ba, bb, kind) or None."""
    fused = _fused_pair(info)
    return None if fused is None else _q_from_fused(info, *fused)

def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])

def _is_methyl(mol: Mol, s: int, ring: set[int]) -> bool:
    return _is_methyl_on_ring(mol, s, ring)

def _methyl_starts_ok(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl(mol, s, ring) for s in starts)

def _q_subs_ok(mol: Mol, ring: set[int], cap: int = 2) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > cap:
        return False
    return True if not starts else _methyl_starts_ok(mol, ring, starts)

def _q_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)

def _is_simple_q_kind(info: dict, want: str) -> bool:
    parts = _q_core(info)
    if parts is None or parts[5] != want or _q_fg_block(info):
        return False
    mol, ring = info["mol"], _ring_set(parts)
    return _outside_ok(mol, ring) and _q_subs_ok(mol, ring)

def _walk_path(ring: list[int], start: int, end: int) -> list[int] | None:
    if start not in ring or end not in ring:
        return None
    i, ln = ring.index(start), len(ring)
    path: list[int] = []
    for k in range(1, ln):
        atom = ring[(i + k) % ln]
        if atom == end:
            return path
        path.append(atom)
    return None

def _path_of_len(ring: list[int], start: int, end: int, n: int) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _walk_path(base, start, end)
        if p is not None and len(p) == n:
            return p
    return None

def _bridge_nb_of(mol: Mol, atom: int, bridge: set[int]) -> list[int]:
    return [
        a.GetIdx() for a in mol.GetAtomWithIdx(atom).GetNeighbors()
        if a.GetIdx() in bridge
    ]

def _quinoline_chain(
    mol: Mol, r1: list[int], r2: list[int], n: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC: 1=N, 2, 3, 4, 4a, 5, 6, 7, 8, 8a (8a = bridge next to N)."""
    nbs = _bridge_nb_of(mol, n, {ba, bb})
    if len(nbs) != 1:
        return None
    a8a, a4a = nbs[0], (bb if nbs[0] == ba else ba)
    pyr, ben = (r1, r2) if n in r1 else (r2, r1)
    mid, ext = _path_of_len(pyr, n, a4a, 3), _path_of_len(ben, a4a, a8a, 4)
    return None if not mid or not ext else [n] + mid + [a4a] + ext + [a8a]

def _iso_a1_a8a(mol: Mol, nbs: list[int], bridge: set[int]) -> tuple[int, int, int] | None:
    """Return (a1, a3, a8a): a1 is N-nb adjacent to bridge a8a."""
    for a1 in nbs:
        hit = _bridge_nb_of(mol, a1, bridge)
        if len(hit) == 1:
            a3 = nbs[0] if a1 == nbs[1] else nbs[1]
            return a1, a3, hit[0]
    return None

def _pyr_nbs(mol: Mol, n: int, pyr: list[int]) -> list[int]:
    pset = set(pyr)
    return [a.GetIdx() for a in mol.GetAtomWithIdx(n).GetNeighbors() if a.GetIdx() in pset]

def _iso_pack(
    mol: Mol, pyr: list[int], ben: list[int], n: int, ba: int, bb: int, nbs: list[int],
) -> list[int] | None:
    got = _iso_a1_a8a(mol, nbs, {ba, bb})
    if got is None:
        return None
    a1, a3, a8a = got
    a4a = bb if a8a == ba else ba
    mid, ext = _path_of_len(pyr, a3, a4a, 1), _path_of_len(ben, a4a, a8a, 4)
    return None if not mid or not ext else [a1, n, a3, mid[0], a4a] + ext + [a8a]

def _isoquinoline_chain(
    mol: Mol, r1: list[int], r2: list[int], n: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC: 1=C next to N, 2=N, 3, 4, 4a, 5..8, 8a."""
    pyr, ben = (r1, r2) if n in r1 else (r2, r1)
    nbs = _pyr_nbs(mol, n, pyr)
    return None if len(nbs) != 2 else _iso_pack(mol, pyr, ben, n, ba, bb, nbs)

def _build_chain(info: dict, parts: tuple) -> list[int]:
    r1, r2, n, ba, bb, kind = parts
    mol = info["mol"]
    if kind == "quinoline":
        return _quinoline_chain(mol, r1, r2, n, ba, bb) or []
    return _isoquinoline_chain(mol, r1, r2, n, ba, bb) or []

def _q_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _q_core(info)
    assert parts is not None
    chain = _build_chain(info, parts)
    return {
        "chain": chain, "n_carbons": 10, "kind": kind,
        "scaffold_id": kind, "n_idx": parts[2], "bridge": [parts[3], parts[4]],
        **extra,
    }

def _try_quinoline_parent(info: dict) -> dict | None:
    if not _is_simple_q_kind(info, "quinoline"):
        return None
    return _q_parent_dict(info, "quinoline")

def _try_isoquinoline_parent(info: dict) -> dict | None:
    if not _is_simple_q_kind(info, "isoquinoline"):
        return None
    return _q_parent_dict(info, "isoquinoline")

def _ol_conflict(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_ester",
        "has_amide", "has_nitrile", "has_amine", "has_thiol",
        "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)

def _is_simple_quinolinol(info: dict) -> bool:
    parts = _q_core(info)
    if parts is None or parts[5] != "quinoline" or _ol_conflict(info):
        return False
    mol, ring = info["mol"], _ring_set(parts)
    oh = _mono_oh_on_ring(info, ring)
    if oh is None or not _outside_ok(mol, ring, {oh["o_idx"]}):
        return False
    return _q_subs_ok(mol, ring)

def _try_quinolinol_parent(info: dict) -> dict | None:
    if not _is_simple_quinolinol(info):
        return None
    return _q_parent_dict(
        info, "quinolinol", oh_c_idx=info["hydroxyls"][0]["c_idx"],
    )

def _is_simple_quinolinecarboxylic(info: dict) -> bool:
    parts = _q_core(info)
    if parts is None or parts[5] != "quinoline":
        return False
    if _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    mol, ring = info["mol"], _ring_set(parts)
    if _carboxyl_ring_c(info, ring) is None:
        return False
    fg_c = info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))

def _try_quinolinecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_quinolinecarboxylic(info):
        return None
    parts = _q_core(info)
    assert parts is not None
    ring = _ring_set(parts)
    return _q_parent_dict(
        info, "quinolinecarboxylic",
        cooh_c_idx=info["carboxyls"][0]["c_idx"],
        ring_attach_idx=_carboxyl_ring_c(info, ring),
    )
