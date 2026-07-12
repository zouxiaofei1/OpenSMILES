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
        if any(n.GetAtomicNum() not in (1, 6) for n in atom.GetNeighbors()):
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
    if am["c_idx"] not in ring_set:
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
