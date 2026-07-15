"""Retained 1H-indazole parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[c]pyrazole. NH=1, N=2 (adjacent);
≤2 methyl/halo; mono CN → carbonitrile; mono CHO → carbaldehyde.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _aldehyde_ring_c,
    _arene_fg_conflict,
    _arene_subs_ok,
    _fg_ring_c,
    _nitrile_n_idx,
)
from namepredict.layer2.fused56 import (
    _all_aromatic,
    _chain_atoms,
    _fused56,
    _six_all_c,
)
from namepredict.layer2.ring_parent import (
    _dbl_o_idx,
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


def _iz_fg_ring(info: dict) -> tuple[Mol, set[int]] | None:
    parts = _iz_parts(info)
    if parts is None:
        return None
    return info["mol"], _ring_set(parts)


def _iz_cn_blocked(info: dict) -> bool:
    allow = frozenset({"has_nitrile"})
    bad = _arene_fg_conflict(info, "has_acid", "has_aldehyde", "has_ketone", allow=allow)
    return not _is_iz_core(info) or bad


def _iz_cn_ctx(info: dict) -> tuple | None:
    if _iz_cn_blocked(info):
        return None
    got = _iz_fg_ring(info)
    if got is None or _fg_ring_c(info, got[1], "nitriles") is None:
        return None
    mol, ring, fg_c = got[0], got[1], info["nitriles"][0]["c_idx"]
    n_idx = _nitrile_n_idx(mol, fg_c)
    return (mol, ring, fg_c, n_idx) if n_idx is not None else None


def _is_simple_indazolecarbonitrile(info: dict) -> bool:
    got = _iz_cn_ctx(info)
    if got is None:
        return False
    mol, ring, fg_c, n_idx = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {n_idx})


def _try_indazolecarbonitrile_parent(info: dict) -> dict | None:
    if not _is_simple_indazolecarbonitrile(info):
        return None
    parts = _iz_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    return _iz_parent_dict(
        info, "indazolecarbonitrile",
        nitrile_c_idx=info["nitriles"][0]["c_idx"],
        ring_attach_idx=_fg_ring_c(info, ring, "nitriles"),
    )


def _iz_ald_ctx(info: dict) -> tuple | None:
    if not _is_iz_core(info) or _arene_fg_conflict(info, "has_acid", "has_ketone"):
        return None
    got = _iz_fg_ring(info)
    if got is None or _aldehyde_ring_c(info, got[1]) is None:
        return None
    mol, ring = got
    fg_c, o_idx = info["aldehydes"][0]["c_idx"], _dbl_o_idx(mol, info["aldehydes"][0]["c_idx"])
    return (mol, ring, fg_c, o_idx) if o_idx is not None else None


def _is_simple_indazolecarbaldehyde(info: dict) -> bool:
    got = _iz_ald_ctx(info)
    if got is None:
        return False
    mol, ring, fg_c, o_idx = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {o_idx})


def _try_indazolecarbaldehyde_parent(info: dict) -> dict | None:
    if not _is_simple_indazolecarbaldehyde(info):
        return None
    parts = _iz_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    return _iz_parent_dict(
        info, "indazolecarbaldehyde",
        aldehyde_c_idx=info["aldehydes"][0]["c_idx"],
        ring_attach_idx=_aldehyde_ring_c(info, ring),
    )
