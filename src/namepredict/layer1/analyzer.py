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


def _is_carboxyl_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    return _has_double_bonded_o(atom) and _has_oh_neighbor(atom)


def _carbon_neighbor_count(atom) -> int:
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])


def _is_ketone_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _carbon_neighbor_count(atom) != 2:
        return False
    return True


def _is_aldehyde_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_double_bonded_o(atom):
        return False
    if _has_oh_neighbor(atom) or _carbon_neighbor_count(atom) > 1:
        return False
    return True


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


def _aldehyde_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if _is_aldehyde_carbon(atom):
            out.append({"c_idx": atom.GetIdx()})
    return out


def _carbon_ids(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


def _fg_flags(
    hydroxyls: list[dict],
    carboxyls: list[dict],
    ketones: list[dict],
    aldehydes: list[dict],
) -> dict:
    return {
        "hydroxyls": hydroxyls,
        "carboxyls": carboxyls,
        "ketones": ketones,
        "aldehydes": aldehydes,
        "has_alcohol": bool(hydroxyls),
        "has_acid": bool(carboxyls),
        "has_ketone": bool(ketones),
        "has_aldehyde": bool(aldehydes),
    }


def _info(mol: Mol, carbons: list[int], fgs: dict) -> dict:
    base = {"mol": mol, "carbon_ids": carbons, "n_carbons": len(carbons)}
    return {**base, **fgs}


def analyze(mol: Mol) -> dict:
    carbons = _carbon_ids(mol)
    fgs = _fg_flags(
        _hydroxyl_entries(mol),
        _carboxyl_entries(mol),
        _ketone_entries(mol),
        _aldehyde_entries(mol),
    )
    return _info(mol, carbons, fgs)
