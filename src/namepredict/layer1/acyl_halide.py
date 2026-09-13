"""L1 检测酰卤 R–C(=O)–X（X = F/Cl/Br/I），依据 IUPAC P-65.5。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C, HALO_Z, O
from namepredict.layer1._carbonyl_common import (
    _alkoxy_c_of,
    _amide_n_of,
    _double_bonded_o_idxs,
    _ester_alkoxy_of as _ester_alkoxy_of_common,
    _has_acid_o_neighbor,
    _has_double_bonded_o,
)

def _is_ester_alkoxy_o(oxygen, carbonyl) -> bool:
    """判断 O 是否为酰卤碳上的中性烷氧基氧。"""
    if oxygen.GetAtomicNum() != O or oxygen.GetTotalNumHs() != 0:
        return False
    if oxygen.GetFormalCharge() != 0:
        return False
    return _alkoxy_c_of(oxygen, carbonyl) is not None

def _ester_alkoxy_of(carbon) -> tuple[int, int] | None:
    """在碳上查找酯样烷氧基侧并返回 (o_idx, alkoxy_c_idx)。"""
    return _ester_alkoxy_of_common(carbon, _is_ester_alkoxy_o)

def _acyl_hal_of(carbon) -> tuple[int, int] | None:
    """返回 F/Cl/Br/I 邻居的 (hal_idx, hal_z)；否则返回 None。"""
    for n in carbon.GetNeighbors():
        z = n.GetAtomicNum()
        if z in HALO_Z:  # 酰卤检测覆盖 F/Cl/Br/I（P-65.5）
            return n.GetIdx(), z
    return None

def _is_acyl_halide_carbon(atom) -> bool:
    """判断碳原子是否为酰卤羰基碳（排除酸、酯、酰胺）。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom) or _ester_alkoxy_of(atom) is not None:
        return False
    if _amide_n_of(atom) is not None:
        return False
    return _acyl_hal_of(atom) is not None

def _entry(atom) -> dict:
    """为酰卤羰基碳组装条目 dict（羰基碳为中心，羰基氧与卤素为周边）。"""
    h = _acyl_hal_of(atom)
    assert h is not None
    hal_idx, _ = h
    return {"center_idx": atom.GetIdx(), "surr_idx": [*_double_bonded_o_idxs(atom), hal_idx]}

def acyl_halide_entries(mol: Mol) -> list[dict]:
    """收集分子中所有酰卤条目的列表。"""
    return [_entry(a) for a in mol.GetAtoms() if _is_acyl_halide_carbon(a)]

def acyl_hal_of(carbon) -> tuple[int, int] | None:
    """返回碳上卤素邻居的 (hal_idx, hal_z)；无则返回 None。"""
    return _acyl_hal_of(carbon)
