from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_alkyl import (
    _disjoint_cover, _is_cf3_carbon, _is_cf3_fluoro, _side_covers, _side_sets,
)


def _all_carbons_are_c(mol: Mol, atom_ids: tuple) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in atom_ids)


def _bond_between(mol: Mol, a: int, b: int):
    return mol.GetBondBetweenAtoms(a, b)


def _ring_bonds_single(mol: Mol, atom_ids: tuple) -> bool:
    ids = list(atom_ids)
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is None or bond.GetBondType().name != "SINGLE":
            return False
    return True


def _outside_carbons(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    return [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
        and a.GetIdx() not in ring_set
        and a.GetIdx() not in skip
    ]


def _pure_c_bonds(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        z = other.GetAtomicNum()
        if z not in (1, 6) or (z != 1 and bond.GetBondType().name != "SINGLE"):
            return False
    return True


def _pure_alkyl_outside(mol: Mol, outside: list[int]) -> bool:
    return all(_is_cf3_carbon(mol, i) or _pure_c_bonds(mol, i) for i in outside)


def _ring_side_starts(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    starts: list[int] = []
    for r in ring_set:
        for n in mol.GetAtomWithIdx(r).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in ring_set:
                if n.GetIdx() not in skip:
                    starts.append(n.GetIdx())
    return starts


def _is_ring_halo(atom, ring_set: set[int]) -> bool:
    if atom.GetAtomicNum() not in (9, 17, 35, 53):
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() in ring_set


def _ring_nitro_n(info: dict, ring_set: set[int]) -> int:
    return sum(1 for n in (info.get("nitros") or []) if n["c_idx"] in ring_set)


def _ring_nitro_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for n in info.get("nitros") or []:
        if n["c_idx"] in ring_set:
            out.add(n["n_idx"])
            out.update(n.get("o_idxs") or [])
    return out


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]

def _outer_fwd(mol: Mol, cur: int, prev: int) -> list:
    atom = mol.GetAtomWithIdx(cur)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return []
    return [x for x in _heavies(atom) if x.GetIdx() != prev]

def _outer_alkoxy_n(mol: Mol, start: int, o_idx: int) -> int:
    fwd = _outer_fwd(mol, start, o_idx)
    if not fwd:
        return 1
    if len(fwd) == 1 and fwd[0].GetAtomicNum() == 6:
        return 2 if not _outer_fwd(mol, fwd[0].GetIdx(), start) else 0
    return 0

def _outer_atoms(mol: Mol, start: int, o_idx: int, n: int) -> list[int]:
    if n == 1:
        return [start]
    fwd = _outer_fwd(mol, start, o_idx)
    return [start, fwd[0].GetIdx()] if fwd else [start]

def _ring_alkoxy_pair(e: dict, ring_set: set[int]) -> tuple[int, int, int] | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in ring_set) == (c2 in ring_set):
        return None
    return (o, c1, c2) if c1 in ring_set else (o, c2, c1)

def _ring_alkoxy_one(mol: Mol, e: dict, ring_set: set[int]) -> dict | None:
    pair = _ring_alkoxy_pair(e, ring_set)
    if pair is None:
        return None
    o, ring_c, outer = pair
    n = _outer_alkoxy_n(mol, outer, o)
    if n == 0:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "n": n,
            "atoms": _outer_atoms(mol, outer, o, n)}

def _ring_alkoxy_ethers(info: dict, ring_set: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _ring_alkoxy_one(mol, e, ring_set)
        if one is not None:
            out.append(one)
    return out

def _ring_alkoxy_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for a in _ring_alkoxy_ethers(info, ring_set):
        out.add(a["o_idx"])
        out.update(a["atoms"])
    return out

def _ring_alkoxy_n(info: dict, ring_set: set[int]) -> int:
    return len(_ring_alkoxy_ethers(info, ring_set))


def _outside_hetero_ok(atom, ring_set: set[int], allow: set[int]) -> bool:
    if atom.GetAtomicNum() in (1, 6) or atom.GetIdx() in ring_set:
        return True
    return atom.GetIdx() in allow or _is_ring_halo(atom, ring_set) or _is_cf3_fluoro(atom)


def _no_hetero_outside(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    allow = allowed or set()
    return all(_outside_hetero_ok(a, ring_set, allow) for a in mol.GetAtoms())


def _outside_ok(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    if not _no_hetero_outside(mol, ring_set, allowed):
        return False
    return _pure_alkyl_outside(mol, _outside_carbons(mol, ring_set, allowed or set()))


def _is_cycloalkane_core(info: dict) -> bool:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return False
    atom_ids = rings[0]["atom_ids"]
    return _all_carbons_are_c(info["mol"], atom_ids) and _ring_bonds_single(info["mol"], atom_ids)


def _is_simple_cycloalkane(info: dict) -> bool:
    if not _is_cycloalkane_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    return _outside_ok(mol, ring_set)


def _ring_double_count(mol: Mol, atom_ids: tuple) -> int:
    ids = list(atom_ids)
    n = 0
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is not None and bond.GetBondType().name == "DOUBLE":
            n += 1
    return n


def _is_carbocycle_ring(info: dict) -> tuple | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    mol: Mol = info["mol"]
    atom_ids = rings[0]["atom_ids"]
    if not _all_carbons_are_c(mol, atom_ids):
        return None
    return atom_ids


def _is_cycloalkene_core(info: dict) -> bool:
    atom_ids = _is_carbocycle_ring(info)
    if atom_ids is None:
        return False
    return _ring_double_count(info["mol"], atom_ids) == 1


def _endocyclic_double(info: dict, ring_set: set[int]) -> tuple[int, int] | None:
    bonds = info.get("double_bonds") or []
    if len(bonds) != 1:
        return None
    c1, c2 = bonds[0]["c1"], bonds[0]["c2"]
    if c1 in ring_set and c2 in ring_set:
        return c1, c2
    return None


def _is_simple_cycloalkene(info: dict) -> bool:
    if not _is_cycloalkene_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if _endocyclic_double(info, ring_set) is None:
        return False
    if not _outside_ok(mol, ring_set):
        return False
    return not _outside_carbons(mol, ring_set)


def _hetero_allowed(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in (1, 6) or atom.GetIdx() in ring_set:
            continue
        if atom.GetIdx() not in allowed:
            return False
    return True


def _mono_oh_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 1:
        return None
    oh = hydroxyls[0]
    if oh["c_idx"] not in ring_set:
        return None
    return oh


def _is_simple_cycloalcohol(info: dict) -> bool:
    if not _is_cycloalkane_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    oh = _mono_oh_on_ring(info, ring_set)
    if oh is None:
        return False
    if not _hetero_allowed(mol, ring_set, {oh["o_idx"]}):
        return False
    return not _outside_carbons(mol, ring_set)


def _mono_amine_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    amines = info.get("amines") or []
    if len(amines) != 1:
        return None
    am = amines[0]
    c_idx = am.get("c_idx")
    if c_idx is None or c_idx not in ring_set:
        return None
    return am


def _is_simple_cycloamine(info: dict) -> bool:
    if not _is_cycloalkane_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    am = _mono_amine_on_ring(info, ring_set)
    if am is None:
        return False
    if not _hetero_allowed(mol, ring_set, {am["n_idx"]}):
        return False
    return not _outside_carbons(mol, ring_set)


def _dbl_o_idx(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != "DOUBLE":
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == 8:
            return other.GetIdx()
    return None


def _mono_ketone_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    ketones = info.get("ketones") or []
    if len(ketones) != 1:
        return None
    ket = ketones[0]
    if ket["c_idx"] not in ring_set:
        return None
    return ket


def _ketone_o_allowed(mol: Mol, ring_set: set[int], ket: dict) -> bool:
    o_idx = _dbl_o_idx(mol, ket["c_idx"])
    if o_idx is None:
        return False
    return _hetero_allowed(mol, ring_set, {o_idx})


def _is_simple_cycloketone(info: dict) -> bool:
    if not _is_cycloalkane_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    ket = _mono_ketone_on_ring(info, ring_set)
    if ket is None or not _ketone_o_allowed(mol, ring_set, ket):
        return False
    return not _outside_carbons(mol, ring_set)


def _is_benzene_core(info: dict) -> bool:
    atom_ids = _is_carbocycle_ring(info)
    if atom_ids is None or len(atom_ids) != 6:
        return False
    mol: Mol = info["mol"]
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)


def _unsub_phenyl_at(mol: Mol, c_idx: int, n_idx: int) -> bool:
    """True if c_idx is the sole N-attachment of an unsubstituted phenyl ring."""
    hits = [set(r) for r in mol.GetRingInfo().AtomRings() if c_idx in r and len(r) == 6]
    if len(hits) != 1: return False
    r = hits[0]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in r): return False
    for i in r:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            z, j = nb.GetAtomicNum(), nb.GetIdx()
            if z != 1 and j not in r and (i != c_idx or j != n_idx or z != 7): return False
    return True


def _ring_halo_n(mol: Mol, ring_set: set[int]) -> int:
    return sum(1 for a in mol.GetAtoms() if _is_ring_halo(a, ring_set))


def _benzene_alkyl_ns(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> list[int]:
    outside = set(_outside_carbons(mol, ring_set, exclude))
    sets = _side_sets(mol, ring_set, starts)
    if sets is None or not _disjoint_cover(sets, outside):
        return []
    return [len(s) for s in sets]


def _mono_benzene_ok(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> bool:
    if not starts:
        return True
    if len(starts) != 1:
        return False
    outside = set(_outside_carbons(mol, ring_set, exclude))
    return _side_covers(mol, starts[0], ring_set, outside)


def _multi_benzene_ok(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> bool:
    ns = _benzene_alkyl_ns(mol, ring_set, starts, exclude)
    return len(ns) == len(starts)


def _benzene_subs_ok(
    mol: Mol, ring_set: set[int], n_nitro: int = 0, n_alkoxy: int = 0,
    exclude: set[int] | None = None,
) -> bool:
    h = _ring_halo_n(mol, ring_set)
    starts = _ring_side_starts(mol, ring_set, exclude)
    n_sub = h + len(starts) + n_nitro + n_alkoxy
    if n_sub > 4:
        return False
    if n_sub <= 1:
        return _mono_benzene_ok(mol, ring_set, starts, exclude)
    return _multi_benzene_ok(mol, ring_set, starts, exclude)


def _arene_alkoxy(info: dict, ring_set: set[int]) -> tuple[set[int], int]:
    alk = _ring_alkoxy_atoms(info, ring_set)
    return alk, len(_ring_alkoxy_ethers(info, ring_set))

def _is_simple_benzene(info: dict) -> bool:
    if not _is_benzene_core(info):
        return False
    mol, ring_set = info["mol"], set(info["rings"][0]["atom_ids"])
    alk, n_alk = _arene_alkoxy(info, ring_set)
    allowed = _ring_nitro_atoms(info, ring_set) | alk
    if not _outside_ok(mol, ring_set, allowed):
        return False
    return _benzene_subs_ok(mol, ring_set, _ring_nitro_n(info, ring_set), n_alk, alk)

def _hetero_or_ring_halo(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    return all(_outside_hetero_ok(a, ring_set, allowed) for a in mol.GetAtoms())

def _arene_fg_subs_ok(
    mol: Mol, ring_set: set[int], allowed: set[int], n_nitro: int = 0,
    n_amino: int = 0, n_alkoxy: int = 0, exclude: set[int] | None = None,
) -> bool:
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    starts = _ring_side_starts(mol, ring_set, exclude)
    if len(_benzene_alkyl_ns(mol, ring_set, starts, exclude)) != len(starts):
        return False
    return _ring_halo_n(mol, ring_set) + len(starts) + n_nitro + n_amino + n_alkoxy <= 4

def _ring_primary_amines(info: dict, ring_set: set[int]) -> list:
    return [
        a for a in (info.get("amines") or [])
        if a.get("degree") == 1 and a.get("c_idx") in ring_set
    ]

def _phenol_amines_ok(info: dict, ring_set: set[int]) -> list | None:
    ring_ams = _ring_primary_amines(info, ring_set)
    if len(info.get("amines") or []) != len(ring_ams):
        return None
    return ring_ams

def _phenol_allowed(info: dict, ring_set: set[int], oh: dict, ring_ams: list) -> set[int]:
    am_n = {a["n_idx"] for a in ring_ams}
    return {oh["o_idx"]} | _ring_nitro_atoms(info, ring_set) | am_n | _ring_alkoxy_atoms(info, ring_set)

def _phenol_subs_ok(mol: Mol, info: dict, ring_set: set[int], oh: dict, ring_ams: list) -> bool:
    alk, n_alk = _arene_alkoxy(info, ring_set)
    allowed = _phenol_allowed(info, ring_set, oh, ring_ams)
    return _arene_fg_subs_ok(
        mol, ring_set, allowed, _ring_nitro_n(info, ring_set), len(ring_ams), n_alk, alk,
    )

def _is_simple_phenol(info: dict) -> bool:
    if not _is_benzene_core(info):
        return False
    mol, ring_set = info["mol"], set(info["rings"][0]["atom_ids"])
    oh, ring_ams = _mono_oh_on_ring(info, ring_set), _phenol_amines_ok(info, ring_set)
    if oh is None or ring_ams is None:
        return False
    return _phenol_subs_ok(mol, info, ring_set, oh, ring_ams)

def _di_oh_on_ring(info: dict, ring_set: set[int]) -> list[dict] | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 2 or any(h["c_idx"] not in ring_set for h in hydroxyls):
        return None
    return hydroxyls

def _is_simple_benzenediol(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("amines"):
        return False
    mol, ring_set = info["mol"], set(info["rings"][0]["atom_ids"])
    ohs = _di_oh_on_ring(info, ring_set)
    if ohs is None or _outside_carbons(mol, ring_set):
        return False
    return _hetero_allowed(mol, ring_set, {h["o_idx"] for h in ohs})

def _benzenediol_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    ohs = _di_oh_on_ring(info, set(ring)) or []
    return {"chain": ring, "n_carbons": len(ring), "kind": "benzenediol",
            "oh_c_idxs": [h["c_idx"] for h in ohs]}

def _is_simple_aniline(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("hydroxyls"):
        return False
    mol, ring_set = info["mol"], set(info["rings"][0]["atom_ids"])
    am = _mono_amine_on_ring(info, ring_set)
    if am is None or am.get("degree") != 1:
        return False
    alk, n_alk = _arene_alkoxy(info, ring_set)
    allowed = {am["n_idx"]} | _ring_nitro_atoms(info, ring_set) | alk
    return _arene_fg_subs_ok(mol, ring_set, allowed, _ring_nitro_n(info, ring_set), 0, n_alk, alk)

def _is_methyl_on_ring(mol: Mol, s: int, ring_set: set[int]) -> bool:
    if _is_cf3_carbon(mol, s):
        return True
    atom = mol.GetAtomWithIdx(s)
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or (z == 6 and n.GetIdx() in ring_set):
            continue
        return False
    return True
