from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import O

def _o_idx(mol: Mol, c_idx: int, bond_name: str) -> int | None:
    """Index of O neighbour of `c_idx` across a bond of the given type."""
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != bond_name:
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == O:
            return other.GetIdx()
    return None


def _dbl_o_idx(mol: Mol, c_idx: int) -> int | None:
    return _o_idx(mol, c_idx, "DOUBLE")
