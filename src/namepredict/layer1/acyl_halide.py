"""L1 检测酰卤 R–C(=O)–X（X = Cl/Br），依据 IUPAC P-65.5。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import Br, C, Cl, O
from namepredict.layer1._carbonyl_common import (
    _alkoxy_c_of,
    _amide_n_of,
    _ester_alkoxy_of as _ester_alkoxy_of_common,
    _has_acid_o_neighbor,
    _has_double_bonded_o,
)

# 酰卤检测仅覆盖 Cl/Br（P-65.5）；F/I 下游不处理
_HAL_Z = frozenset({Cl, Br})

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
    """返回 Cl/Br 邻居的 (hal_idx, hal_z)；否则返回 None。"""
    for n in carbon.GetNeighbors():
        z = n.GetAtomicNum()
        if z in _HAL_Z:
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
    """为酰卤羰基碳组装条目 dict（含卤索引与原子序数）。"""
    h = _acyl_hal_of(atom)
    assert h is not None
    hal_idx, hal_z = h
    return {
        "c_idx": atom.GetIdx(),
        "hal_idx": hal_idx,
        "hal_z": hal_z,
        "cl_idx": hal_idx,  # 兼容：L2/L3 过滤器仍使用 cl_idx
    }

def acyl_halide_entries(mol: Mol) -> list[dict]:
    """收集分子中所有酰卤条目的列表。"""
    return [_entry(a) for a in mol.GetAtoms() if _is_acyl_halide_carbon(a)]

def acyl_hal_of(carbon) -> tuple[int, int] | None:
    """返回碳上卤素邻居的 (hal_idx, hal_z)；无则返回 None。"""
    return _acyl_hal_of(carbon)
