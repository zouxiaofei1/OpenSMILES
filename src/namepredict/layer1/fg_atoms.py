"""L1 官能团特征原子"""
from __future__ import annotations

from namepredict.constants import O


def _idx(payload: dict, key: str) -> int | None:
    """取 payload 中单个索引键（缺失返回 None）。"""
    v = payload.get(key)
    return int(v) if v is not None else None


def _hetero_neighbors(mol, idx: int | None, z: int) -> set[int]:
    """返回指定原子的某种元素邻居索引。"""
    if idx is None:
        return set()
    return {n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetAtomicNum() == z}


def _phosphate_atoms(mol, payload: dict) -> set[int]:
    """磷酸：P 中心 + 其全部氧（=O 与三个单键 O，含 O–R 桥氧）。"""
    p = _idx(payload, "p_idx")
    return (set() if p is None else {p}) | _hetero_neighbors(mol, p, O)


def center_surr_atoms(payload: dict) -> frozenset[int]:
    """通用 FG 特征原子：中心原子与全部周边原子之并。"""
    out = {int(i) for i in (payload.get("surr_idx") or ())}
    center = _idx(payload, "center_idx")
    if center is not None:
        out.add(center)
    return frozenset(out)


FG_ATOM_FNS = {  # 不走通用规则的例外类别
    "phosphate": _phosphate_atoms,
}
