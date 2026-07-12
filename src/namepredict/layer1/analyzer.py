from __future__ import annotations

from rdkit.Chem import Mol


def _is_hydroxyl_oxygen(atom) -> bool:
    if atom.GetAtomicNum() != 8:
        return False
    if atom.GetTotalNumHs() < 1:
        return False
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(carbons) == 1


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


def _carbon_ids(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


def _info(mol: Mol, carbons: list[int], hydroxyls: list[dict]) -> dict:
    return {
        "mol": mol,
        "carbon_ids": carbons,
        "n_carbons": len(carbons),
        "hydroxyls": hydroxyls,
        "has_alcohol": bool(hydroxyls),
    }


def analyze(mol: Mol) -> dict:
    carbons = _carbon_ids(mol)
    hydroxyls = _hydroxyl_entries(mol)
    return _info(mol, carbons, hydroxyls)
