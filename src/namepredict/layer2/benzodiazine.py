"""Retained quinazoline parent (IUPAC P-22.2.1 / P-25).

Fused 6+6 aromatic 8C+2N (pyrimidine fusion: N-1 adjacent to bridge, N-3 meta).
Unsubstituted only (this round).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.naphthalene import (
    _bridge_adjacent,
    _bridge_pair,
    _two_six_rings,
)
from namepredict.layer2.ring_parent import _outside_ok, _ring_halo_n, _ring_side_starts


def _fused_pair(info: dict):
    pair = _two_six_rings(info)
    if pair is None:
        return None
    br = _bridge_pair(*pair)
    if br is None or not _bridge_adjacent(info["mol"], *br):
        return None
    return (*pair, br)


def _atoms_ok(mol: Mol, r1: list[int], r2: list[int]) -> set[int] | None:
    atoms = set(r1) | set(r2)
    if len(atoms) != 10:
        return None
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms):
        return None
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atoms]
    return atoms if zs.count(7) == 2 and zs.count(6) == 8 else None


def _n_idxs(mol: Mol, atoms: set[int]) -> list[int]:
    return sorted(i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)


def _n_adj_bridge(mol: Mol, n: int, bridge: set[int]) -> bool:
    nbs = {a.GetIdx() for a in mol.GetAtomWithIdx(n).GetNeighbors()}
    return bool(nbs & bridge)


def _n_share_mid(mol: Mol, n0: int, n1: int, bridge: set[int]) -> bool:
    mid = (
        {a.GetIdx() for a in mol.GetAtomWithIdx(n0).GetNeighbors()}
        & {a.GetIdx() for a in mol.GetAtomWithIdx(n1).GetNeighbors()}
    )
    return len(mid - bridge) >= 1


def _is_quinazoline_n(mol: Mol, ns: list[int], bridge: set[int]) -> bool:
    if any(n in bridge for n in ns):
        return False
    adj = sum(1 for n in ns if _n_adj_bridge(mol, n, bridge))
    return adj == 1 and _n_share_mid(mol, ns[0], ns[1], bridge)


def _qz_from_fused(mol: Mol, r1, r2, br) -> tuple[set[int], list[int], set[int]] | None:
    atoms = _atoms_ok(mol, r1, r2)
    if atoms is None:
        return None
    ns = _n_idxs(mol, atoms)
    if not _is_quinazoline_n(mol, ns, set(br)):
        return None
    return atoms, ns, set(br)


def _qz_core(info: dict) -> tuple[set[int], list[int], set[int]] | None:
    fused = _fused_pair(info)
    return None if fused is None else _qz_from_fused(info["mol"], *fused)


def _unsub_ok(mol: Mol, atoms: set[int]) -> bool:
    return (
        _outside_ok(mol, atoms)
        and _ring_halo_n(mol, atoms) == 0
        and not _ring_side_starts(mol, atoms)
    )


def _fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_ether", "has_acyl_chloride",
        "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _is_simple_quinazoline(info: dict) -> bool:
    core = _qz_core(info)
    if core is None or _fg_block(info):
        return False
    return _unsub_ok(info["mol"], core[0])


def _n1(mol: Mol, ns: list[int], bridge: set[int]) -> int:
    for n in ns:
        nbs = {a.GetIdx() for a in mol.GetAtomWithIdx(n).GetNeighbors()}
        if nbs & bridge:
            return n
    return ns[0]


def _quinazoline_parent(info: dict) -> dict:
    core = _qz_core(info)
    atoms, ns, bridge = core
    n1 = _n1(info["mol"], ns, bridge)
    chain = [n1] + [i for i in sorted(atoms) if i != n1]
    return {
        "chain": chain, "n_carbons": 10, "kind": "quinazoline",
        "ring_atoms": list(atoms),
    }


def _try_quinazoline_parent(info: dict) -> dict | None:
    return _quinazoline_parent(info) if _is_simple_quinazoline(info) else None
