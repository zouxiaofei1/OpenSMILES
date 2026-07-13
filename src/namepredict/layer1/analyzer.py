from __future__ import annotations

from rdkit.Chem import BondType, Mol

def _is_single_c_oh(atom) -> bool:
    if atom.GetAtomicNum() != 8 or atom.GetTotalNumHs() < 1:
        return False
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) == 1

def _dbl_o_on(bond, carbon) -> bool:
    if bond.GetBondType() != BondType.DOUBLE:
        return False
    return bond.GetOtherAtom(carbon).GetAtomicNum() == 8

def _has_double_bonded_o(carbon) -> bool:
    return any(_dbl_o_on(b, carbon) for b in carbon.GetBonds())

def _has_oh_neighbor(carbon) -> bool:
    return any(_is_single_c_oh(n) for n in carbon.GetNeighbors())

def _acyl_cl_of(carbon) -> int | None:
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() == 17:
            return n.GetIdx()
    return None

def _is_acyl_chloride_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _ester_alkoxy_of(atom) is not None:
        return False
    if _amide_n_of(atom) is not None:
        return False
    return _acyl_cl_of(atom) is not None

def _is_carboxyl_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    return _has_double_bonded_o(atom) and _has_oh_neighbor(atom)

def _carbon_neighbor_count(atom) -> int:
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])

def _amide_n_info(carbon) -> tuple[int, list[int]] | None:
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != 7:
            continue
        o = [x for x in n.GetNeighbors()
             if x.GetAtomicNum() != 1 and x.GetIdx() != carbon.GetIdx()]
        if len(o) <= 2 and all(x.GetAtomicNum() == 6 for x in o):
            return n.GetIdx(), [x.GetIdx() for x in o]
    return None

def _amide_n_of(carbon) -> int | None:
    info = _amide_n_info(carbon)
    return info[0] if info else None

def _is_amide_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _ester_alkoxy_of(atom) is not None:
        return False
    return _amide_n_info(atom) is not None

def _is_ketone_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _carbon_neighbor_count(atom) != 2:
        return False
    if _amide_n_of(atom) is not None:
        return False
    if _anhydride_o_of(atom) is not None:
        return False
    return True

def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == 6 and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None

def _is_anhydride_bridge_o(oxygen) -> bool:
    if oxygen.GetAtomicNum() != 8 or oxygen.GetTotalNumHs() != 0:
        return False
    cs = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(cs) != 2:
        return False
    return all(_has_double_bonded_o(c) and not _has_oh_neighbor(c) for c in cs)

def _anhydride_o_of(carbon) -> int | None:
    for n in carbon.GetNeighbors():
        if _is_anhydride_bridge_o(n):
            return n.GetIdx()
    return None

def _ether_cs(oxygen) -> list:
    return [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() == 6]

def _is_ether_oxygen(atom) -> bool:
    if atom.GetAtomicNum() != 8 or atom.GetTotalNumHs() != 0:
        return False
    if _is_anhydride_bridge_o(atom):
        return False
    cs = _ether_cs(atom)
    return len(cs) == 2 and not any(_has_double_bonded_o(c) for c in cs)

def _ether_entry(atom) -> dict:
    cs = _ether_cs(atom)
    c1, c2 = cs[0].GetIdx(), cs[1].GetIdx()
    return {"o_idx": atom.GetIdx(), "c1": c1, "c2": c2}

def _ether_entries(mol: Mol) -> list[dict]:
    return [_ether_entry(a) for a in mol.GetAtoms() if _is_ether_oxygen(a)]

def _sulfide_cs(sulfur) -> list:
    return [n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6]

def _is_sulfide_sulfur(atom) -> bool:
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 2:
        return False
    cs = _sulfide_cs(atom)
    return len(cs) == 2 and not any(_has_double_bonded_o(c) for c in cs)

def _sulfide_entry(atom) -> dict:
    cs = _sulfide_cs(atom)
    return {"s_idx": atom.GetIdx(), "c1": cs[0].GetIdx(), "c2": cs[1].GetIdx()}

def _sulfide_entries(mol: Mol) -> list[dict]:
    return [_sulfide_entry(a) for a in mol.GetAtoms() if _is_sulfide_sulfur(a)]

def _is_nitro_nitrogen(atom) -> bool:
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
        return False
    nbs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    os_ = [n for n in nbs if n.GetAtomicNum() == 8]
    cs = [n for n in nbs if n.GetAtomicNum() == 6]
    return len(nbs) == 3 and len(os_) == 2 and len(cs) == 1

def _nitro_entry(atom) -> dict:
    nbs = list(atom.GetNeighbors())
    o_idxs = [n.GetIdx() for n in nbs if n.GetAtomicNum() == 8]
    c_idx = next(n.GetIdx() for n in nbs if n.GetAtomicNum() == 6)
    return {"n_idx": atom.GetIdx(), "c_idx": c_idx, "o_idxs": o_idxs}

def _nitro_entries(mol: Mol) -> list[dict]:
    return [_nitro_entry(a) for a in mol.GetAtoms() if _is_nitro_nitrogen(a)]

def _is_ester_alkoxy_o(oxygen, carbonyl) -> bool:
    if oxygen.GetAtomicNum() != 8 or oxygen.GetTotalNumHs() != 0:
        return False
    if _is_anhydride_bridge_o(oxygen):
        return False
    return _alkoxy_c_of(oxygen, carbonyl) is not None

def _ester_alkoxy_of(carbon) -> tuple[int, int] | None:
    for n in carbon.GetNeighbors():
        if not _is_ester_alkoxy_o(n, carbon):
            continue
        alkoxy = _alkoxy_c_of(n, carbon)
        if alkoxy is not None:
            return n.GetIdx(), alkoxy
    return None

def _is_ester_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom):
        return False
    return _ester_alkoxy_of(atom) is not None

def _is_aldehyde_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _carbon_neighbor_count(atom) > 1:
        return False
    if _ester_alkoxy_of(atom) is not None or _acyl_cl_of(atom) is not None:
        return False
    if _anhydride_o_of(atom) is not None:
        return False
    return _amide_n_of(atom) is None

def _is_hydroxyl_oxygen(atom) -> bool:
    if not _is_single_c_oh(atom):
        return False
    return not _is_carboxyl_carbon(_carbon_neighbor(atom))

def _carbon_neighbor(atom):
    return next(n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)

def _hydroxyl_entry(atom) -> dict:
    carbon = _carbon_neighbor(atom)
    return {"o_idx": atom.GetIdx(), "c_idx": carbon.GetIdx()}

def _hydroxyl_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_hydroxyl_oxygen(atom):
            out.append(_hydroxyl_entry(atom))
    return out

def _thiol_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() < 1:
            continue
        if _carbon_neighbor_count(atom) != 1:
            continue
        out.append({"s_idx": atom.GetIdx(), "c_idx": _carbon_neighbor(atom).GetIdx()})
    return out

def _is_amide_n(atom) -> bool:
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() == 6 and _has_double_bonded_o(n):
            return True
    return False

def _amine_degree(atom) -> int | None:
    if atom.GetAtomicNum() != 7 or _is_amide_n(atom):
        return None
    n_c, n_h = _carbon_neighbor_count(atom), atom.GetTotalNumHs()
    if n_c == 1 and n_h >= 2:
        return 1
    return 2 if n_c == 2 and n_h == 1 else (3 if n_c == 3 and n_h == 0 else None)

def _amine_entry(atom, deg: int) -> dict:
    cs = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    base = {"n_idx": atom.GetIdx(), "degree": deg}
    return {**base, "c_idxs": cs} if deg >= 2 else {**base, "c_idx": cs[0]}

def _amine_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        deg = _amine_degree(atom)
        if deg is not None:
            out.append(_amine_entry(atom, deg))
    return out

def _carboxyl_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_carboxyl_carbon(atom):
            out.append({"c_idx": atom.GetIdx()})
    return out

def _ketone_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_ketone_carbon(atom):
            out.append({"c_idx": atom.GetIdx()})
    return out

def _amide_entry(atom) -> dict:
    n_idx, n_cs = _amide_n_info(atom)
    return {"c_idx": atom.GetIdx(), "n_idx": n_idx, "n_c_idxs": n_cs}

def _amide_entries(mol: Mol) -> list[dict]:
    return [
        _amide_entry(a) for a in mol.GetAtoms() if _is_amide_carbon(a)
    ]

def _aldehyde_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_aldehyde_carbon(atom):
            out.append({"c_idx": atom.GetIdx()})
    return out

def _acyl_chloride_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_acyl_chloride_carbon(atom):
            out.append({"c_idx": atom.GetIdx(), "cl_idx": _acyl_cl_of(atom)})
    return out

def _ester_entry(atom) -> dict:
    o_idx, alkoxy_c = _ester_alkoxy_of(atom)
    return {"c_idx": atom.GetIdx(), "o_idx": o_idx, "alkoxy_c_idx": alkoxy_c}

def _ester_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_ester_carbon(atom):
            out.append(_ester_entry(atom))
    return out

def _anhydride_other_c(oxygen, carbon) -> int:
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == 6 and n.GetIdx() != carbon.GetIdx():
            return n.GetIdx()
    return carbon.GetIdx()

def _anhydride_entry(atom) -> dict:
    o_idx = _anhydride_o_of(atom)
    o_atom = atom.GetOwningMol().GetAtomWithIdx(o_idx)
    other = _anhydride_other_c(o_atom, atom)
    c1, c2 = sorted((atom.GetIdx(), other))
    return {"o_idx": o_idx, "c1_idx": c1, "c2_idx": c2}

def _is_anhydride_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom):
        return False
    return _anhydride_o_of(atom) is not None

def _anhydride_key(e: dict) -> tuple[int, int, int]:
    return e["o_idx"], e["c1_idx"], e["c2_idx"]

def _anhydride_entries(mol: Mol) -> list[dict]:
    seen: set[tuple[int, int, int]] = set()
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if not _is_anhydride_carbon(atom):
            continue
        e = _anhydride_entry(atom)
        if _anhydride_key(e) not in seen:
            seen.add(_anhydride_key(e))
            out.append(e)
    return out

def _is_cc_double(bond) -> bool:
    if bond.GetBondType() != BondType.DOUBLE or bond.GetIsAromatic():
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6

def _is_cc_triple(bond) -> bool:
    if bond.GetBondType() != BondType.TRIPLE:
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6

def _bond_entry(bond) -> dict:
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    return {"c1": min(a, b), "c2": max(a, b)}

def _double_bond_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for bond in mol.GetBonds():
        if _is_cc_double(bond):
            out.append(_bond_entry(bond))
    return out

def _triple_bond_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for bond in mol.GetBonds():
        if _is_cc_triple(bond):
            out.append(_bond_entry(bond))
    return out

def _is_cn_triple(bond) -> bool:
    if bond.GetBondType() != BondType.TRIPLE:
        return False
    z = {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
    return z == {6, 7}

def _nitrile_entry(bond) -> dict:
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    c = a if a.GetAtomicNum() == 6 else b
    n = b if a.GetAtomicNum() == 6 else a
    return {"c_idx": c.GetIdx(), "n_idx": n.GetIdx()}

def _nitrile_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for bond in mol.GetBonds():
        if _is_cn_triple(bond):
            out.append(_nitrile_entry(bond))
    return out

def _ring_entry(atom_ids: tuple) -> dict:
    return {"atom_ids": atom_ids}

def _ring_entries(mol: Mol) -> list[dict]:
    return [_ring_entry(r) for r in mol.GetRingInfo().AtomRings()]

def _ring_meta(mol: Mol) -> dict:
    rings = _ring_entries(mol)
    return {"rings": rings, "n_rings": len(rings), "has_ring": bool(rings)}

def _carbon_ids(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]

def _fg_core_lists(
    hydroxyls: list[dict],
    carboxyls: list[dict],
    esters: list[dict],
    amides: list[dict],
    ketones: list[dict],
) -> dict:
    return {
        "hydroxyls": hydroxyls,
        "carboxyls": carboxyls,
        "esters": esters,
        "amides": amides,
        "ketones": ketones,
    }

def _fg_more_lists(parts: dict) -> dict:
    keys = (
        "aldehydes", "amines", "nitriles", "double_bonds", "triple_bonds",
        "acyl_chlorides", "anhydrides", "thiols", "ethers", "sulfides",
        "nitros",
    )
    return {k: parts[k] for k in keys}

def _fg_lists(parts: dict) -> dict:
    core = _fg_core_lists(
        parts["hydroxyls"], parts["carboxyls"], parts["esters"],
        parts["amides"], parts["ketones"],
    )
    return {**core, **_fg_more_lists(parts)}

def _fg_bools_core(lists: dict) -> dict:
    return {
        "has_alcohol": bool(lists["hydroxyls"]),
        "has_acid": bool(lists["carboxyls"]),
        "has_ester": bool(lists["esters"]),
        "has_amide": bool(lists["amides"]),
        "has_ketone": bool(lists["ketones"]),
    }

def _fg_bools_more(lists: dict) -> dict:
    keys = (
        ("has_aldehyde", "aldehydes"), ("has_amine", "amines"),
        ("has_nitrile", "nitriles"), ("has_alkene", "double_bonds"),
        ("has_alkyne", "triple_bonds"), ("has_acyl_chloride", "acyl_chlorides"),
        ("has_anhydride", "anhydrides"), ("has_thiol", "thiols"),
        ("has_ether", "ethers"), ("has_sulfide", "sulfides"),
        ("has_nitro", "nitros"),
    )
    return {hk: bool(lists[lk]) for hk, lk in keys}

def _fg_bools(lists: dict) -> dict:
    return {**_fg_bools_core(lists), **_fg_bools_more(lists)}

def _fg_parts_a(mol: Mol) -> dict:
    return {
        "hydroxyls": _hydroxyl_entries(mol),
        "carboxyls": _carboxyl_entries(mol),
        "esters": _ester_entries(mol),
        "amides": _amide_entries(mol),
        "ketones": _ketone_entries(mol),
    }

def _fg_parts_b(mol: Mol) -> dict:
    return {
        "aldehydes": _aldehyde_entries(mol), "amines": _amine_entries(mol),
        "nitriles": _nitrile_entries(mol), "double_bonds": _double_bond_entries(mol),
        "triple_bonds": _triple_bond_entries(mol),
        "acyl_chlorides": _acyl_chloride_entries(mol),
        "anhydrides": _anhydride_entries(mol), "thiols": _thiol_entries(mol),
        "ethers": _ether_entries(mol), "sulfides": _sulfide_entries(mol),
        "nitros": _nitro_entries(mol),
    }

def _fg_parts(mol: Mol) -> dict:
    return {**_fg_parts_a(mol), **_fg_parts_b(mol)}

def _collect_fgs(mol: Mol) -> dict:
    lists = _fg_lists(_fg_parts(mol))
    return {**lists, **_fg_bools(lists)}

def _info(mol: Mol, carbons: list[int], fgs: dict) -> dict:
    base = {"mol": mol, "carbon_ids": carbons, "n_carbons": len(carbons)}
    return {**base, **fgs, **_ring_meta(mol)}

def analyze(mol: Mol) -> dict:
    carbons = _carbon_ids(mol)
    return _info(mol, carbons, _collect_fgs(mol))
