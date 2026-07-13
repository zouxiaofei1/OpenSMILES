from __future__ import annotations

from rdkit.Chem import Mol


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


def _outside_carbons(mol: Mol, ring_set: set[int]) -> list[int]:
    return [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring_set
    ]


def _pure_alkyl_outside(mol: Mol, outside: list[int]) -> bool:
    for idx in outside:
        atom = mol.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            z = other.GetAtomicNum()
            if z not in (1, 6) or (z != 1 and bond.GetBondType().name != "SINGLE"):
                return False
    return True


def _ring_side_starts(mol: Mol, ring_set: set[int]) -> list[int]:
    starts: list[int] = []
    for r in ring_set:
        for n in mol.GetAtomWithIdx(r).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in ring_set:
                starts.append(n.GetIdx())
    return starts


def _is_ring_halo(atom, ring_set: set[int]) -> bool:
    if atom.GetAtomicNum() not in (9, 17, 35, 53):
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() in ring_set


def _no_hetero_outside(mol: Mol, ring_set: set[int]) -> bool:
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in (1, 6) or atom.GetIdx() in ring_set:
            continue
        if _is_ring_halo(atom, ring_set):
            continue
        return False
    return True


def _outside_ok(mol: Mol, ring_set: set[int]) -> bool:
    if not _no_hetero_outside(mol, ring_set):
        return False
    return _pure_alkyl_outside(mol, _outside_carbons(mol, ring_set))


def _is_cycloalkane_core(info: dict) -> bool:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return False
    mol: Mol = info["mol"]
    atom_ids = rings[0]["atom_ids"]
    if not _all_carbons_are_c(mol, atom_ids):
        return False
    return _ring_bonds_single(mol, atom_ids)


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


def _ring_halo_n(mol: Mol, ring_set: set[int]) -> int:
    return sum(1 for a in mol.GetAtoms() if _is_ring_halo(a, ring_set))


def _benzene_alkyl_ns(mol: Mol, ring_set: set[int], starts: list[int]) -> list[int]:
    outside = set(_outside_carbons(mol, ring_set))
    if len(outside) != len(starts):
        return []
    if any(s not in outside for s in starts):
        return []
    return [1] * len(starts)


def _mono_benzene_ok(mol: Mol, ring_set: set[int], starts: list[int]) -> bool:
    if not starts:
        return True
    return len(_outside_carbons(mol, ring_set)) in (1, 2)


def _multi_benzene_ok(mol: Mol, ring_set: set[int], starts: list[int]) -> bool:
    ns = _benzene_alkyl_ns(mol, ring_set, starts)
    return len(ns) == len(starts) and all(n == 1 for n in ns)


def _benzene_subs_ok(mol: Mol, ring_set: set[int]) -> bool:
    h = _ring_halo_n(mol, ring_set)
    starts = _ring_side_starts(mol, ring_set)
    n_sub = h + len(starts)
    if n_sub > 3:
        return False
    if n_sub <= 1:
        return _mono_benzene_ok(mol, ring_set, starts)
    return _multi_benzene_ok(mol, ring_set, starts)


def _is_simple_benzene(info: dict) -> bool:
    if not _is_benzene_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring_set):
        return False
    return _benzene_subs_ok(mol, ring_set)


def _hetero_or_ring_halo(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in (1, 6) or atom.GetIdx() in ring_set:
            continue
        if atom.GetIdx() in allowed or _is_ring_halo(atom, ring_set):
            continue
        return False
    return True


def _arene_fg_subs_ok(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    starts = _ring_side_starts(mol, ring_set)
    if len(_benzene_alkyl_ns(mol, ring_set, starts)) != len(starts):
        return False
    return _ring_halo_n(mol, ring_set) + len(starts) <= 2


def _is_simple_phenol(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("amines"):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    oh = _mono_oh_on_ring(info, ring_set)
    if oh is None:
        return False
    return _arene_fg_subs_ok(mol, ring_set, {oh["o_idx"]})


def _di_oh_on_ring(info: dict, ring_set: set[int]) -> list[dict] | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 2:
        return None
    if any(h["c_idx"] not in ring_set for h in hydroxyls):
        return None
    return hydroxyls


def _is_simple_benzenediol(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("amines"):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    ohs = _di_oh_on_ring(info, ring_set)
    if ohs is None or _outside_carbons(mol, ring_set):
        return False
    return _hetero_allowed(mol, ring_set, {h["o_idx"] for h in ohs})


def _benzenediol_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    ohs = _di_oh_on_ring(info, set(ring)) or []
    return {
        "chain": ring,
        "n_carbons": len(ring),
        "kind": "benzenediol",
        "oh_c_idxs": [h["c_idx"] for h in ohs],
    }


def _is_simple_aniline(info: dict) -> bool:
    if not _is_benzene_core(info) or info.get("hydroxyls"):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    am = _mono_amine_on_ring(info, ring_set)
    if am is None or am.get("degree") != 1:
        return False
    return _arene_fg_subs_ok(mol, ring_set, {am["n_idx"]})


def _ring_c_neighbors(mol: Mol, c_idx: int, ring_set: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() in ring_set
    ]


def _carboxyl_ring_c(info: dict, ring_set: set[int]) -> int | None:
    carboxyls = info.get("carboxyls") or []
    if len(carboxyls) != 1:
        return None
    nbs = _ring_c_neighbors(info["mol"], carboxyls[0]["c_idx"], ring_set)
    return nbs[0] if len(nbs) == 1 else None


def _cooh_oxygen_idxs(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8}


def _ring_phenol_ohs(info: dict, ring_set: set[int]) -> list[dict]:
    return [h for h in (info.get("hydroxyls") or []) if h["c_idx"] in ring_set]


def _benzoic_conflict_fg(info: dict) -> bool:
    bad = (
        "has_ester", "has_amide", "has_aldehyde", "has_ketone",
        "has_acyl_chloride", "has_anhydride", "has_nitrile",
    )
    if any(info.get(k) for k in bad):
        return True
    return bool(info.get("amines") or info.get("thiols") or info.get("ethers"))


def _is_methyl_on_ring(mol: Mol, s: int, ring_set: set[int]) -> bool:
    atom = mol.GetAtomWithIdx(s)
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or (z == 6 and n.GetIdx() in ring_set):
            continue
        return False
    return True


def _benzoic_alkyl_ok(mol: Mol, ring_set: set[int], cooh_c: int) -> bool:
    starts = [s for s in _ring_side_starts(mol, ring_set) if s != cooh_c]
    outside = [i for i in _outside_carbons(mol, ring_set) if i != cooh_c]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring_set) for s in starts)


def _benzoic_extra_n(info: dict, mol: Mol, ring_set: set[int], cooh_c: int) -> int:
    starts = [s for s in _ring_side_starts(mol, ring_set) if s != cooh_c]
    ohs = _ring_phenol_ohs(info, ring_set)
    return _ring_halo_n(mol, ring_set) + len(starts) + len(ohs)


def _benzoic_hetero_ok(info: dict, mol: Mol, ring_set: set[int], cooh_c: int) -> bool:
    ohs = _ring_phenol_ohs(info, ring_set)
    allowed = _cooh_oxygen_idxs(mol, cooh_c) | {h["o_idx"] for h in ohs}
    return _hetero_or_ring_halo(mol, ring_set, allowed)


def _benzoic_subs_ok(info: dict, mol: Mol, ring_set: set[int], cooh_c: int) -> bool:
    if not _benzoic_hetero_ok(info, mol, ring_set, cooh_c):
        return False
    if not _benzoic_alkyl_ok(mol, ring_set, cooh_c):
        return False
    return _benzoic_extra_n(info, mol, ring_set, cooh_c) <= 2


def _is_simple_benzoic(info: dict) -> bool:
    if not _is_benzene_core(info) or _benzoic_conflict_fg(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if _carboxyl_ring_c(info, ring_set) is None:
        return False
    return _benzoic_subs_ok(info, mol, ring_set, info["carboxyls"][0]["c_idx"])


def _benzoic_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 6, "kind": "benzoic",
        "cooh_c_idx": cooh_c, "ring_attach_idx": _carboxyl_ring_c(info, set(ring)),
    }


def _try_benzoic_parent(info: dict) -> dict | None:
    return _benzoic_parent(info) if _is_simple_benzoic(info) else None
