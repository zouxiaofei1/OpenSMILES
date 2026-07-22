"""L1 detection of aryl/alkyl boronic acids R–B(OH)2 (P-68.1)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import B, C, H, O


def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]


def _bond_type(a, b):
    bond = a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
    return bond.GetBondType() if bond is not None else None


def _is_sgl(a, b) -> bool:
    return _bond_type(a, b) == BondType.SINGLE


def _is_oh_on_b(oxygen, boron) -> bool:
    """Neutral O–H single-bonded only to B (not ether/oxide)."""
    if oxygen.GetAtomicNum() != O or oxygen.GetFormalCharge() != 0:
        return False
    if oxygen.GetTotalNumHs() < 1 or oxygen.GetTotalDegree() != 2:
        return False
    return _is_sgl(oxygen, boron)


def _oh_nbs(boron) -> list:
    return [n for n in _heavies(boron) if _is_oh_on_b(n, boron)]


def _c_nbs(boron) -> list:
    return [n for n in boron.GetNeighbors() if n.GetAtomicNum() == C]


def _is_boronic_b(atom) -> bool:
    """Tricoordinate neutral B with two OH and one C."""
    if atom.GetAtomicNum() != B or atom.GetFormalCharge() != 0:
        return False
    if atom.GetTotalDegree() != 3 or len(_heavies(atom)) != 3:
        return False
    if len(_oh_nbs(atom)) != 2:
        return False
    return len(_c_nbs(atom)) == 1


def _boronic_entry(atom) -> dict:
    ohs = _oh_nbs(atom)
    c = _c_nbs(atom)[0]
    return {
        "b_idx": atom.GetIdx(),
        "c_attach": c.GetIdx(),
        "o_idxs": [o.GetIdx() for o in ohs],
    }


def boronic_entries(mol: Mol) -> list[dict]:
    return [_boronic_entry(a) for a in mol.GetAtoms() if _is_boronic_b(a)]
