"""末端母体所有权：不可变 owned_atoms 的最终确定。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _dbl_o_idx, _o_idx


def _chain_atoms(parent: dict) -> set[int]:
    """取母体链原子集合。"""
    return set(parent.get("chain") or [])


def _add_opt(out: set[int], idx: int | None) -> set[int]:
    """非空索引加入集合并返回。"""
    if idx is not None:
        out.add(idx)
    return out


def _amide_n_from_c(mol: Mol, c_idx: int) -> int | None:
    """找连接碳的酰胺氮索引。"""
    carbon = mol.GetAtomWithIdx(c_idx)
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _amide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """酰胺：羰基 C + =O + N。"""
    c_idx = parent.get("amide_c_idx")
    if c_idx is None:
        return set()
    out = {c_idx}
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _amide_n_from_c(mol, c_idx))


def _single_o_idx(mol: Mol, c_idx: int) -> int | None:
    """跨单键找碳的 O 邻居索引。"""
    return _o_idx(mol, c_idx, "SINGLE")


def _acid_o_atoms(mol: Mol, c_idx: int) -> set[int]:
    """羧酸碳的 =O 与 -O- 氧原子。"""
    out: set[int] = set()
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _single_o_idx(mol, c_idx))


def _cooh_c_idxs(parent: dict) -> list[int]:
    """归一化取母体的 cooh 碳索引列表。"""
    multi = parent.get("cooh_c_idxs")
    if multi:
        return [int(x) for x in multi]
    one = parent.get("cooh_c_idx")
    return [int(one)] if one is not None else []


def _acid_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """羧酸各 cooh 碳的 FG 原子并集。"""
    out: set[int] = set()
    for c in _cooh_c_idxs(parent):
        out.add(c)
        out |= _acid_o_atoms(mol, c)
    return out


def _aldehyde_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """醛：羰基 C + =O。"""
    c_idx = parent.get("aldehyde_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    return _add_opt(out, _dbl_o_idx(mol, int(c_idx)))


def _ether_arm_atoms(mol: Mol, o_idx: int, c_idx: int) -> set[int]:
    """醚氧一侧的碳臂原子集（禁走 O）。"""
    from namepredict.layer2.chain_walk import _longest_from
    return {c_idx, *_longest_from(mol, c_idx, {o_idx})}


def _chalcogen_arm_fg(mol: Mol, parent: dict, kind: str, idx_field: str) -> set[int]:
    """硫族中心 + 碳臂（短臂不一定在链中）。"""
    if parent.get("kind") != kind or parent.get(idx_field) is None:
        return set()
    center = int(parent[idx_field])
    out = {center}
    for n in mol.GetAtomWithIdx(center).GetNeighbors():
        if n.GetAtomicNum() == 6:
            out |= _ether_arm_atoms(mol, center, n.GetIdx())
    return out


def _ether_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """醚 O + 两条碳臂（短臂不一定在链中）。"""
    return _chalcogen_arm_fg(mol, parent, "ether", "o_idx")


def _sulfide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """硫醚 S + 两条碳臂（短臂不一定在链中）。"""
    return _chalcogen_arm_fg(mol, parent, "sulfide", "s_idx")


def _hydroxy_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """醇/二醇/三醇：连接碳 + OH 氧。"""
    c_idxs = parent.get("oh_c_idxs") or ([parent.get("oh_c_idx")] if parent.get("oh_c_idx") is not None else [])
    if not c_idxs:
        return set()
    out: set[int] = set()
    for c in c_idxs:
        out.add(int(c))
        _add_opt(out, _single_o_idx(mol, int(c)))
    return out


def _ketone_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """酮羰基 C + 双键 O（多酮取全部羰基；存在时含乙酰甲基）。"""
    c_idxs = parent.get("ketone_c_idxs") or ([parent.get("ketone_c_idx")] if parent.get("ketone_c_idx") is not None else [])
    if not c_idxs:
        return set()
    out: set[int] = set()
    for c in c_idxs:
        out.add(int(c))
        _add_opt(out, _dbl_o_idx(mol, int(c)))
    return _add_opt(out, parent.get("acetyl_methyl_idx"))


def _amine_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """苯胺/胺：母体链连接碳 + 胺氮（N 上非母体臂留在所有权外作 N- 取代基）。"""
    c_idxs = parent.get("amine_c_idxs") or ([parent.get("amine_c_idx")] if parent.get("amine_c_idx") is not None else [])
    if not c_idxs:
        return set()
    chain = set(parent.get("chain") or [])
    out: set[int] = set()
    for c in c_idxs:
        if c not in chain:
            continue  # 非母体链的 N 臂：所有权外，作为 N- 取代基
        out.add(int(c))
        for n in mol.GetAtomWithIdx(int(c)).GetNeighbors():
            if n.GetAtomicNum() == 7:
                out.add(n.GetIdx())
    return out


def _nitrile_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """腈/苯甲腈：CN 碳 + 三键氮。"""
    c_idx = parent.get("nitrile_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 7:
            out.add(n.GetIdx())
    return out


def _thiol_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """硫醇/二硫醇：连接碳 + SH 硫。"""
    c_idxs = parent.get("sh_c_idxs") or ([parent.get("sh_c_idx")] if parent.get("sh_c_idx") is not None else [])
    if not c_idxs:
        return set()
    out: set[int] = set()
    for c in c_idxs:
        out.add(int(c))
        for n in mol.GetAtomWithIdx(int(c)).GetNeighbors():
            if n.GetAtomicNum() == 16:
                out.add(n.GetIdx())
    return out


def _one_ester_fg(mol: Mol, c_idx: int, out: set[int]) -> None:
    """添加一个酯羰基 C + =O + -O-（仅酸侧，无烷氧基臂）。"""
    out.add(int(c_idx))
    dbl_o = _dbl_o_idx(mol, int(c_idx))
    _add_opt(out, dbl_o)
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 8 and n.GetIdx() != dbl_o:
            out.add(n.GetIdx())


def _single_ester_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """单酯/苯甲酸酯：羰基 C + =O + -O- + 烷氧基臂。"""
    c_idx = parent.get("ester_c_idx")
    if c_idx is None:
        return set()
    out: set[int] = set()
    _one_ester_fg(mol, c_idx, out)
    return out


def _ester_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """单酯 FG 原子（复用单酯逻辑）。"""
    return _single_ester_fg_atoms(mol, parent) 


def _anhydride_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """酸酐：两个酰基碳 + 桥接 O + 两个羰基氧。"""
    out: set[int] = set()
    for key in ("acyl_c_idx", "other_acyl_c_idx"):
        c = parent.get(key)
        if c is None:
            continue
        out.add(int(c))
        _add_opt(out, _dbl_o_idx(mol, int(c)))
    return _add_opt(out, parent.get("o_idx"))


def _acyl_chloride_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """酰氯：酰基碳 + 羰基 O + Cl。"""
    c_idx = parent.get("acyl_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    _add_opt(out, _dbl_o_idx(mol, int(c_idx)))
    return _add_opt(out, parent.get("cl_idx") if parent.get("cl_idx") is not None else parent.get("hal_idx"))


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    """FG 所有权由字段驱动：拥有主官能团的每个重原子。"""
    parts = (
        _amide_fg_atoms(mol, parent) if parent.get("amide_c_idx") is not None else set(),
        _aldehyde_fg_atoms(mol, parent) if parent.get("aldehyde_c_idx") is not None else set(),
        _acid_fg_atoms(mol, parent),
        _ether_fg_atoms(mol, parent),
        _sulfide_fg_atoms(mol, parent),
        _hydroxy_fg_atoms(mol, parent),
        _ketone_fg_atoms(mol, parent),
        _amine_fg_atoms(mol, parent),
        _nitrile_fg_atoms(mol, parent),
        _thiol_fg_atoms(mol, parent),
        _ester_fg_atoms(mol, parent),
        _anhydride_fg_atoms(mol, parent),
        _acyl_chloride_fg_atoms(mol, parent),
    )
    out: set[int] = set()
    for part in parts:
        out |= part
    return out


def compute_owned_atoms(parent: dict, mol: Mol) -> frozenset[int]:
    """链与 kind 特异的 FG 原子的并集（末端所有权集合）。"""
    return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))


def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
    """一次性复制候选，生成不可变 owned_atoms frozenset。"""
    if isinstance(parent.get("owned_atoms"), frozenset):
        return parent
    return {**parent, "owned_atoms": compute_owned_atoms(parent, mol)}
