"""9,10-Anthraquinone retained parent (IUPAC P-25 / P-64).

Linear fused three C6 rings (14 C) with exactly two meso ring ketones at 9,10.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.anthracene import _is_linear
from namepredict.layer2.scaffold.ring_parent import (
    _dbl_o_idx,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _aq_shape_ok(s: dict) -> bool:
    return (
        s.get("n_rings") == 3
        and s.get("n_atoms") == 14
        and s.get("topology") == "fused"
        and not s.get("hetero_atoms")
        and len(s.get("fusion_edges") or []) == 2
    )


def _aq_system(info: dict) -> dict | None:
    for s in info.get("ring_systems") or []:
        if _aq_shape_ok(s) and _is_linear(s):
            return s
    return None


def _ketone_idxs(info: dict) -> list[int]:
    return [e["c_idx"] for e in (info.get("ketones") or []) if "c_idx" in e]


def _center_and_outers(info: dict) -> tuple[list[int], list[list[int]]] | None:
    ket = _ketone_idxs(info)
    if len(ket) != 2:
        return None
    rings = [list(r["atom_ids"]) for r in (info.get("rings") or [])]
    center = next((r for r in rings if set(ket) <= set(r)), None)
    if center is None:
        return None
    outers = [r for r in rings if set(r) & set(center) and r is not center]
    return (center, outers) if len(outers) == 2 else None


def _c_nbs(mol: Mol, idx: int) -> list[int]:
    return [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetAtomicNum() == 6]


def _is_bridgehead(info: dict, n: int) -> bool:
    return sum(1 for r in info["rings"] if n in r["atom_ids"]) == 2


def _ketone_meso_one(info: dict, mol: Mol, k: int) -> bool:
    nbs = _c_nbs(mol, k)
    return len(nbs) == 2 and all(_is_bridgehead(info, n) for n in nbs)


def _ketones_are_meso(info: dict, mol: Mol) -> bool:
    """Both ketones on the shared center ring; each C-nb is a fusion bridgehead."""
    if _center_and_outers(info) is None:
        return False
    return all(_ketone_meso_one(info, mol, k) for k in _ketone_idxs(info))


def _ketone_o_idxs(mol: Mol, ket: list[int]) -> set[int] | None:
    out: set[int] = set()
    for c in ket:
        o = _dbl_o_idx(mol, c)
        if o is None:
            return None
        out.add(o)
    return out


def _is_methyl_c(mol: Mol, s: int) -> bool:
    return mol.GetAtomWithIdx(s).GetDegree() == 1


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    return set(outside) == set(starts) and all(_is_methyl_c(mol, s) for s in starts)


def _aq_subs_ok(mol: Mol, ring: set[int]) -> bool:
    if _ring_halo_n(mol, ring) != 0:
        return False
    starts = _ring_side_starts(mol, ring)
    return not starts or _mono_methyl_only(mol, ring, starts)


def _fg_extra_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_alcohol", "has_ester", "has_amide",
        "has_nitrile", "has_amine", "has_thiol", "has_nitro", "has_ether",
        "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _aq_env_ok(info: dict, mol: Mol, system: dict) -> bool:
    ring = set(system["atom_ids"])
    oset = _ketone_o_idxs(mol, _ketone_idxs(info))
    return oset is not None and _outside_ok(mol, ring, oset) and _aq_subs_ok(mol, ring)


def _is_simple_anthraquinone(info: dict) -> bool:
    system = _aq_system(info)
    if system is None or _fg_extra_block(info) or len(_ketone_idxs(info)) != 2:
        return False
    mol: Mol = info["mol"]
    return _ketones_are_meso(info, mol) and _aq_env_ok(info, mol, system)


def _path_len4(ring: list[int], ba: int, bb: int) -> list[int] | None:
    if ba not in ring:
        return None
    i, path, n = ring.index(ba), [], len(ring)
    for k in range(1, n):
        atom = ring[(i + k) % n]
        if atom == bb:
            return path if len(path) == 4 else None
        path.append(atom)
    return None


def _exterior4(ring: list[int], ba: int, bb: int) -> list[int] | None:
    return _path_len4(ring, ba, bb) or _path_len4(list(reversed(ring)), ba, bb)


def _chain_for_pair(
    mol: Mol, outers: list[list[int]], k9: int, k10: int, n9: list[int], n10: list[int],
) -> list[int] | None:
    a9, a8, a4, a10a = n9[0], n9[1], n10[0], n10[1]
    o14 = next((o for o in outers if a9 in o and a4 in o), None)
    o58 = next((o for o in outers if a8 in o and a10a in o), None)
    if o14 is None or o58 is None or o14 is o58:
        return None
    e14, e58 = _exterior4(o14, a9, a4), _exterior4(o58, a10a, a8)
    if e14 is None or e58 is None:
        return None
    return e14 + [a4, k10, a10a] + e58 + [a8, k9, a9]


def _flip_pairs(nbs: list[int]) -> list[list[int]]:
    return [nbs, list(reversed(nbs))]


def _chains_for_ket_pair(mol: Mol, outers: list[list[int]], k9: int, k10: int) -> list[list[int]]:
    n9, n10 = _c_nbs(mol, k9), _c_nbs(mol, k10)
    if len(n9) != 2 or len(n10) != 2:
        return []
    out: list[list[int]] = []
    for a, b in ((p, q) for p in _flip_pairs(n9) for q in _flip_pairs(n10)):
        c = _chain_for_pair(mol, outers, k9, k10, a, b)
        if c is not None:
            out.append(c)
    return out


def _uniq_chains(raw: list[list[int]]) -> list[list[int]]:
    seen: set[tuple] = set()
    uniq: list[list[int]] = []
    for c in raw:
        t = tuple(c)
        if t not in seen:
            seen.add(t)
            uniq.append(c)
    return uniq


def _aq_chains(info: dict) -> list[list[int]]:
    parts = _center_and_outers(info)
    if parts is None:
        return []
    _, outers = parts
    mol, ket = info["mol"], _ketone_idxs(info)
    raw = _chains_for_ket_pair(mol, outers, ket[0], ket[1])
    raw += _chains_for_ket_pair(mol, outers, ket[1], ket[0])
    return _uniq_chains(raw)


def _anthraquinone_parent(info: dict) -> dict:
    system = _aq_system(info) or {}
    chains = _aq_chains(info)
    chain = chains[0] if chains else list(system.get("atom_ids") or [])
    ket = _ketone_idxs(info)
    return {
        "chain": chain, "n_carbons": 14, "kind": "anthraquinone",
        "scaffold_id": "anthraquinone", "ring_atoms": list(system.get("atom_ids") or []),
        "ketone_c_idxs": ket, "anthra_chains": chains,
    }


def _try_anthraquinone_parent(info: dict) -> dict | None:
    return _anthraquinone_parent(info) if _is_simple_anthraquinone(info) else None
