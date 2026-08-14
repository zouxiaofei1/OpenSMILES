from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import O

def _o_idx(mol: Mol, c_idx: int, bond_name: str) -> int | None:
    """跨给定类型的键，返回 `c_idx` 的 O 邻居索引。"""
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
