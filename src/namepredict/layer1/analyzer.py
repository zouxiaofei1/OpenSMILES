from __future__ import annotations

from rdkit.Chem import Mol


def _is_hydroxyl_oxygen(atom) -> bool:
    if atom.GetAtomicNum() != 8:
        return False
    if atom.GetTotalNumHs() < 1:
        return False
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(carbons) == 1


def _hydroxyl_entries(mol: Mol) -> list[dict]:
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if not _is_hydroxyl_oxygen(atom):
            continue
        carbon = next(n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        out.append({"o_idx": atom.GetIdx(), "c_idx": carbon.GetIdx()})
    return out


def _carbon_ids(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


def analyze(mol: Mol) -> dict:
    carbons = _carbon_ids(mol)
    hydroxyls = _hydroxyl_entries(mol)
    return {
        "mol": mol,
        "carbon_ids": carbons,
        "n_carbons": len(carbons),
        "hydroxyls": hydroxyls,
        "has_alcohol": bool(hydroxyls),
    }
