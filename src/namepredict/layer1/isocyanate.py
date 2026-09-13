"""L1 检测异氰酸酯 R–N=C=O 和异硫氰酸酯 R–N=C=S（P-61.9）。"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, H, N, O, S


def _bond_type(a, b):
    """返回两原子之间的键类型。"""
    bond = a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
    return bond.GetBondType() if bond is not None else None


def _is_dbl(a, b) -> bool:
    """判断两原子之间是否为双键。"""
    return _bond_type(a, b) == BondType.DOUBLE


def _is_sgl(a, b) -> bool:
    """判断两原子之间是否为单键。"""
    return _bond_type(a, b) == BondType.SINGLE


def _heavies(atom):
    """返回原子连有的非氢重原子邻居。"""
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]


def _pick_z(nbs, z: int):
    """从邻居中挑出第一个指定原子序数的原子。"""
    return next((a for a in nbs if a.GetAtomicNum() == z), None)


def _cumul_pair(carbon, x_z: int):
    """若碳为二配位双键结构 N=C=X 则返回 (N, X)；否则返回 None。"""
    if carbon.GetAtomicNum() != C or carbon.GetTotalDegree() != 2:
        return None
    nbs = _heavies(carbon)
    if len(nbs) != 2:
        return None
    n, x = _pick_z(nbs, N), _pick_z(nbs, x_z)
    if n is None or x is None or not (_is_dbl(carbon, n) and _is_dbl(carbon, x)):
        return None
    return n, x


def _r_of_iso_n(n_atom, carbon) -> int | None:
    """与 N 单键相连的 R 碳（而非累积双键的 C）。"""
    for n in _heavies(n_atom):
        if n.GetIdx() == carbon.GetIdx():
            continue
        if n.GetAtomicNum() == C and _is_sgl(n_atom, n):
            return n.GetIdx()
    return None


def _iso_n_ok(n_atom, carbon) -> bool:
    """判断 N 是否为二配位且连有 R 碳的异氰酸酯氮。"""
    if n_atom.GetAtomicNum() != N or n_atom.GetTotalDegree() != 2:
        return False
    return _r_of_iso_n(n_atom, carbon) is not None


def _pack_entry(carbon, n_atom, x) -> dict:
    """组装异氰酸酯条目 dict（累积碳为中心，氮与 X 为周边）。"""
    return {"center_idx": carbon.GetIdx(), "surr_idx": [n_atom.GetIdx(), x.GetIdx()]}


def _entry_for(carbon, x_z: int) -> dict | None:
    """为碳生成异氰酸酯/异硫氰酸酯条目；不匹配返回 None。"""
    pair = _cumul_pair(carbon, x_z)
    if pair is None:
        return None
    n_atom, x = pair
    if not _iso_n_ok(n_atom, carbon):
        return None
    r = _r_of_iso_n(n_atom, carbon)
    return None if r is None else _pack_entry(carbon, n_atom, x)


def _entries(mol: Mol, x_z: int) -> list[dict]:
    """收集分子中指定 X 的异氰酸酯类条目列表。"""
    return [e for a in mol.GetAtoms() if (e := _entry_for(a, x_z)) is not None]


def isocyanate_entries(mol: Mol) -> list[dict]:
    """收集分子中所有异氰酸酯条目（X = O）。"""
    return _entries(mol, O)


def isothiocyanate_entries(mol: Mol) -> list[dict]:
    """收集分子中所有异硫氰酸酯条目（X = S）。"""
    return _entries(mol, S)
