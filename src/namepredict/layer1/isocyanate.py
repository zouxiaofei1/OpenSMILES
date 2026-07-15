"""L1 detection of isocyanate R–N=C=O and isothiocyanate R–N=C=S (P-61.9)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol


def _bond_type(a, b):
    bond = a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
    return bond.GetBondType() if bond is not None else None


def _is_dbl(a, b) -> bool:
    return _bond_type(a, b) == BondType.DOUBLE


def _is_sgl(a, b) -> bool:
    return _bond_type(a, b) == BondType.SINGLE


def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _pick_z(nbs, z: int):
    return next((a for a in nbs if a.GetAtomicNum() == z), None)


def _cumul_pair(carbon, x_z: int):
    """(N, X) if carbon is digonal N=C=X; else None."""
    if carbon.GetAtomicNum() != 6 or carbon.GetTotalDegree() != 2:
        return None
    nbs = _heavies(carbon)
    if len(nbs) != 2:
        return None
    n, x = _pick_z(nbs, 7), _pick_z(nbs, x_z)
    if n is None or x is None or not (_is_dbl(carbon, n) and _is_dbl(carbon, x)):
        return None
    return n, x


def _r_of_iso_n(n_atom, carbon) -> int | None:
    """R carbon single-bonded to N (not the cumulated C)."""
    for n in _heavies(n_atom):
        if n.GetIdx() == carbon.GetIdx():
            continue
        if n.GetAtomicNum() == 6 and _is_sgl(n_atom, n):
            return n.GetIdx()
    return None


def _iso_n_ok(n_atom, carbon) -> bool:
    if n_atom.GetAtomicNum() != 7 or n_atom.GetTotalDegree() != 2:
        return False
    return _r_of_iso_n(n_atom, carbon) is not None


def _pack_entry(carbon, n_atom, x, r: int) -> dict:
    return {
        "c_idx": carbon.GetIdx(), "n_idx": n_atom.GetIdx(),
        "x_idx": x.GetIdx(), "r_c_idx": r,
    }


def _entry_for(carbon, x_z: int) -> dict | None:
    pair = _cumul_pair(carbon, x_z)
    if pair is None:
        return None
    n_atom, x = pair
    if not _iso_n_ok(n_atom, carbon):
        return None
    r = _r_of_iso_n(n_atom, carbon)
    return None if r is None else _pack_entry(carbon, n_atom, x, r)


def _entries(mol: Mol, x_z: int) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry_for(a, x_z)) is not None]


def isocyanate_entries(mol: Mol) -> list[dict]:
    return _entries(mol, 8)


def isothiocyanate_entries(mol: Mol) -> list[dict]:
    return _entries(mol, 16)
