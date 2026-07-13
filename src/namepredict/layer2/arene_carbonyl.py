"""Arene retained carbonyl parents: benzoic acid, benzaldehyde, acetophenone."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _dbl_o_idx,
    _hetero_or_ring_halo,
    _is_benzene_core,
    _is_methyl_on_ring,
    _outside_carbons,
    _ring_halo_n,
    _ring_side_starts,
)


def _ring_c_neighbors(mol: Mol, c_idx: int, ring_set: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() in ring_set
    ]


def _fg_ring_c(info: dict, ring_set: set[int], ekey: str) -> int | None:
    entries = info.get(ekey) or []
    if len(entries) != 1:
        return None
    nbs = _ring_c_neighbors(info["mol"], entries[0]["c_idx"], ring_set)
    return nbs[0] if len(nbs) == 1 else None


def _carboxyl_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "carboxyls")


def _aldehyde_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "aldehydes")


def _ketone_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "ketones")


def _cooh_oxygen_idxs(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8}


def _ring_phenol_ohs(info: dict, ring_set: set[int]) -> list[dict]:
    return [h for h in (info.get("hydroxyls") or []) if h["c_idx"] in ring_set]


def _arene_fg_conflict(info: dict, *extra: str) -> bool:
    bad = (
        "has_ester", "has_amide", "has_acyl_chloride", "has_anhydride",
        "has_nitrile", *extra,
    )
    if any(info.get(k) for k in bad):
        return True
    return bool(info.get("amines") or info.get("thiols") or info.get("ethers"))


def _arene_alkyl_ok(mol: Mol, ring_set: set[int], exclude: set[int]) -> bool:
    starts = [s for s in _ring_side_starts(mol, ring_set) if s not in exclude]
    outside = [i for i in _outside_carbons(mol, ring_set) if i not in exclude]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring_set) for s in starts)


def _arene_extra_n(info: dict, mol: Mol, ring_set: set[int], exclude: set[int]) -> int:
    starts = [s for s in _ring_side_starts(mol, ring_set) if s not in exclude]
    ohs = _ring_phenol_ohs(info, ring_set)
    return _ring_halo_n(mol, ring_set) + len(starts) + len(ohs)


def _arene_hetero_ok(info: dict, mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    ohs = _ring_phenol_ohs(info, ring_set)
    return _hetero_or_ring_halo(mol, ring_set, allowed | {h["o_idx"] for h in ohs})


def _arene_subs_ok(
    info: dict, mol: Mol, ring_set: set[int], exclude: set[int], allowed: set[int],
) -> bool:
    if not _arene_hetero_ok(info, mol, ring_set, allowed):
        return False
    if not _arene_alkyl_ok(mol, ring_set, exclude):
        return False
    return _arene_extra_n(info, mol, ring_set, exclude) <= 2


def _is_simple_benzoic(info: dict) -> bool:
    if not _is_benzene_core(info) or _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if _carboxyl_ring_c(info, ring_set) is None:
        return False
    fg_c = info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring_set, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))


def _benzoic_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "benzoic",
        "cooh_c_idx": cooh_c, "ring_attach_idx": _carboxyl_ring_c(info, set(ring)),
    }


def _try_benzoic_parent(info: dict) -> dict | None:
    return _benzoic_parent(info) if _is_simple_benzoic(info) else None


def _oxo_fg_ok(info: dict, ekey: str, *conflict: str) -> tuple | None:
    if not _is_benzene_core(info) or _arene_fg_conflict(info, *conflict):
        return None
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if _fg_ring_c(info, ring, ekey) is None:
        return None
    fg_c = info[ekey][0]["c_idx"]
    o_idx = _dbl_o_idx(mol, fg_c)
    return (mol, ring, fg_c, o_idx) if o_idx is not None else None


def _is_simple_benzaldehyde(info: dict) -> bool:
    got = _oxo_fg_ok(info, "aldehydes", "has_acid", "has_ketone")
    if got is None:
        return False
    mol, ring, fg_c, o_idx = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {o_idx})


def _benzaldehyde_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    ald_c = info["aldehydes"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "benzaldehyde",
        "aldehyde_c_idx": ald_c, "ring_attach_idx": _aldehyde_ring_c(info, set(ring)),
    }


def _try_benzaldehyde_parent(info: dict) -> dict | None:
    return _benzaldehyde_parent(info) if _is_simple_benzaldehyde(info) else None


def _is_methyl_carbon(mol: Mol, c_idx: int, only_nb: int) -> bool:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == only_nb


def _acetyl_methyl_c(mol: Mol, ket_c: int, ring_set: set[int]) -> int | None:
    atom = mol.GetAtomWithIdx(ket_c)
    cands = [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() not in ring_set
    ]
    if len(cands) != 1:
        return None
    me = cands[0]
    return me if _is_methyl_carbon(mol, me, ket_c) else None


def _is_simple_acetophenone(info: dict) -> bool:
    got = _oxo_fg_ok(info, "ketones", "has_acid", "has_aldehyde")
    if got is None:
        return False
    mol, ring, ket_c, o_idx = got
    me = _acetyl_methyl_c(mol, ket_c, ring)
    if me is None:
        return False
    return _arene_subs_ok(info, mol, ring, {ket_c, me}, {o_idx})


def _acetophenone_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    ket_c = info["ketones"][0]["c_idx"]
    me = _acetyl_methyl_c(info["mol"], ket_c, set(ring))
    return {
        "chain": ring, "n_carbons": 6, "kind": "acetophenone",
        "ketone_c_idx": ket_c, "acetyl_methyl_idx": me,
        "ring_attach_idx": _ketone_ring_c(info, set(ring)),
    }


def _try_acetophenone_parent(info: dict) -> dict | None:
    return _acetophenone_parent(info) if _is_simple_acetophenone(info) else None
