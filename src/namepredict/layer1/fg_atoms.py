"""L1 官能团特征原子：各 FG 类别「由哪些重原子构成」的唯一定义（P-14/P-44 所有权与 P-14.4 位次的共同底料）。

此处只含 FG 的内在原子（由分子结构唯一决定）；「哪些原子进了母体骨架、附着点是谁」属关系信息，
由 L2 结合所选骨架另行推导。下游一律消费 ``occurrence.characteristic_atoms``，不再各自回读 payload。
"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, O


def _idx(payload: dict, key: str) -> int | None:
    """取 payload 中单个索引键（缺失返回 None）。"""
    v = payload.get(key)
    return int(v) if v is not None else None


def _pair(payload: dict, *keys: str) -> set[int]:
    """按多个单值索引键收集非空索引。"""
    return {x for k in keys if (x := _idx(payload, k)) is not None}


def _os_of(mol: Mol, carbon: int | None, bond_type) -> set[int]:
    """返回碳上指定键型的全部氧邻居索引。"""
    if carbon is None:
        return set()
    out: set[int] = set()
    for n in mol.GetAtomWithIdx(carbon).GetNeighbors():
        if n.GetAtomicNum() != O:
            continue
        b = mol.GetBondBetweenAtoms(carbon, n.GetIdx())
        if b is not None and b.GetBondType() == bond_type:
            out.add(n.GetIdx())
    return out


def _carbonyl_pair(mol: Mol, carbon: int | None) -> set[int]:
    """羰基碳与其双键氧（碳缺失时返回空集）。"""
    return set() if carbon is None else {carbon} | _os_of(mol, carbon, BondType.DOUBLE)


def _hetero_neighbors(mol: Mol, idx: int | None, z: int) -> set[int]:
    """返回指定原子的某种元素邻居索引。"""
    if idx is None:
        return set()
    return {n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetAtomicNum() == z}


def _radical_atoms(mol: Mol, payload: dict) -> set[int]:
    """自由基连接点碳（`*` 锚点为虚拟原子，不入母体原子集）。"""
    return _pair(payload, "c_idx")


def _acyl_atoms(mol: Mol, payload: dict) -> set[int]:
    """酰基残基：羰基头碳 + 羰基氧（P-65.1.7.2）。"""
    return _carbonyl_pair(mol, _idx(payload, "c_idx"))


def _acid_atoms(mol: Mol, payload: dict) -> set[int]:
    """羧酸：羧基碳 + 羰基氧 + 羟基氧（或羧酸根氧）。"""
    c = _idx(payload, "c_idx")
    return _carbonyl_pair(mol, c) | _os_of(mol, c, BondType.SINGLE)


def _phosphate_atoms(mol: Mol, payload: dict) -> set[int]:
    """磷酸：P 中心 + 其全部氧（=O 与三个单键 O，含 O–R 桥氧）。"""
    p = _idx(payload, "p_idx")
    return (set() if p is None else {p}) | _hetero_neighbors(mol, p, O)


def _anhydride_atoms(mol: Mol, payload: dict) -> set[int]:
    """酸酐：两个酰基碳 + 各自羰基氧 + 桥氧。"""
    out: set[int] = set()
    for key in ("c1_idx", "c2_idx"):
        out |= _carbonyl_pair(mol, _idx(payload, key))
    return out | _pair(payload, "o_idx")


def _ester_atoms(mol: Mol, payload: dict) -> set[int]:
    """酯：羰基碳 + 羰基氧 + 酯氧（烷氧基臂属取代基侧，不属 FG 本体）。"""
    c = _idx(payload, "c_idx")
    return _carbonyl_pair(mol, c) | _os_of(mol, c, BondType.SINGLE)


def _acyl_halide_atoms(mol: Mol, payload: dict) -> set[int]:
    """酰卤：羰基碳 + 羰基氧 + 卤素。"""
    return _carbonyl_pair(mol, _idx(payload, "c_idx")) | _pair(payload, "hal_idx")


def _amide_atoms(mol: Mol, payload: dict) -> set[int]:
    """酰胺：羰基碳 + 羰基氧 + 酰胺氮（N 上取代基由 L3 作 N- 前缀）。"""
    return _carbonyl_pair(mol, _idx(payload, "c_idx")) | _pair(payload, "n_idx")


def _nitrile_atoms(mol: Mol, payload: dict) -> set[int]:
    """腈：腈碳 + 三键氮。"""
    return _pair(payload, "c_idx", "n_idx")


def _aldehyde_atoms(mol: Mol, payload: dict) -> set[int]:
    """醛：羰基碳 + 羰基氧。"""
    return _carbonyl_pair(mol, _idx(payload, "c_idx"))


def _ketone_atoms(mol: Mol, payload: dict) -> set[int]:
    """酮：羰基碳 + 羰基氧（乙酰甲基属骨架关系，由 L2 判定）。"""
    return _carbonyl_pair(mol, _idx(payload, "c_idx"))


def _alcohol_atoms(mol: Mol, payload: dict) -> set[int]:
    """醇：连接碳 + 羟基氧。"""
    return _pair(payload, "c_idx", "o_idx")


def _thiol_atoms(mol: Mol, payload: dict) -> set[int]:
    """硫醇：连接碳 + 巯基硫。"""
    return _pair(payload, "c_idx", "s_idx")


def _amine_atoms(mol: Mol, payload: dict) -> set[int]:
    """胺：氮 + 全部碳臂（哪些臂在母体链上属骨架关系，由 L2 判定）。"""
    n = _idx(payload, "n_idx")
    return (set() if n is None else {n}) | _hetero_neighbors(mol, n, C)


def _ether_atoms(mol: Mol, payload: dict) -> set[int]:
    """醚：氧 + 两条碳臂。"""
    return _pair(payload, "o_idx", "c1", "c2")


def _sulfide_atoms(mol: Mol, payload: dict) -> set[int]:
    """硫醚：硫 + 两条碳臂。"""
    return _pair(payload, "s_idx", "c1", "c2")


def _quaternary_ammonium_atoms(mol: Mol, payload: dict) -> set[int]:
    """季铵：四配位氮 + 其全部碳。"""
    n = _idx(payload, "n_idx")
    return (set() if n is None else {n}) | _hetero_neighbors(mol, n, C)


def _isocyanate_atoms(mol: Mol, payload: dict) -> set[int]:
    """异氰酸酯/异硫氰酸酯：R–N=C=X 的 N、蓄积碳与 X（O 或 S）。"""
    return _pair(payload, "c_idx", "n_idx", "x_idx")


FG_ATOM_FNS = {  # FG 类别值（FgSpec.fg）→ 特征原子函数；新增 FG 时在此登记，未登记则退回 payload 索引键猜测。
    "radical": _radical_atoms,
    "acyl": _acyl_atoms,
    "acid": _acid_atoms,
    "phosphate": _phosphate_atoms,
    "anhydride": _anhydride_atoms,
    "ester": _ester_atoms,
    "acyl_halide": _acyl_halide_atoms,
    "amide": _amide_atoms,
    "nitrile": _nitrile_atoms,
    "aldehyde": _aldehyde_atoms,
    "ketone": _ketone_atoms,
    "alcohol": _alcohol_atoms,
    "thiol": _thiol_atoms,
    "amine": _amine_atoms,
    "ether": _ether_atoms,
    "sulfide": _sulfide_atoms,
    "quaternary_ammonium": _quaternary_ammonium_atoms,
    "isocyanate": _isocyanate_atoms,
    "isothiocyanate": _isocyanate_atoms,
}
