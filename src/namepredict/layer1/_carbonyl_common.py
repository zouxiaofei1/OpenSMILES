"""layer1 检测器共享的羰基检测原语；检测器各自保留 `_is_ester_alkoxy_o` 谓词并传入共享的 `_ester_alkoxy_of`。
"""
from __future__ import annotations

from rdkit.Chem import BondType

from namepredict.constants import C, H, N, O


def _is_single_c_oh(atom) -> bool:
    """判断 O 原子是否为与单一碳相连的羟基氧。"""
    if atom.GetAtomicNum() != O or atom.GetTotalNumHs() < 1:
        return False
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == C]) == 1


def _dbl_o_on(bond, carbon) -> bool:
    """判断 bond 是否为指向 O 的双键。"""
    if bond.GetBondType() != BondType.DOUBLE:
        return False
    return bond.GetOtherAtom(carbon).GetAtomicNum() == O


def _has_double_bonded_o(carbon) -> bool:
    """判断碳是否连有羰基双键氧。"""
    return any(_dbl_o_on(b, carbon) for b in carbon.GetBonds())


def _double_bonded_o_idxs(carbon) -> list[int]:
    """返回碳上羰基双键氧的索引列表。"""
    mol = carbon.GetOwningMol()
    out: list[int] = []
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != O:
            continue
        b = mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
        if b is not None and b.GetBondType() == BondType.DOUBLE:
            out.append(n.GetIdx())
    return out


def _is_carboxylate_o(atom) -> bool:
    """判断 O 是否为羧酸盐阴离子氧。"""
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != -1:
        return False
    return atom.GetTotalDegree() == 1 and atom.GetTotalNumHs() == 0


def _has_acid_o_neighbor(carbon) -> bool:
    """判断碳是否连有酸性羟基或羧酸盐氧邻居。"""
    return any(
        _is_single_c_oh(n) or _is_carboxylate_o(n) for n in carbon.GetNeighbors()
    )


def _is_anhydride_bridge_o(oxygen) -> bool:
    """判断 O 是否为连接两个羰基的酸酐桥氧。"""
    if oxygen.GetAtomicNum() != O or oxygen.GetTotalNumHs() != 0:
        return False
    cs = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() == C]
    if len(cs) != 2:
        return False
    return all(_has_double_bonded_o(c) and not _has_acid_o_neighbor(c) for c in cs)


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    """返回氧上除羰基碳外的烷氧基碳索引。"""
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


def _ester_alkoxy_of(carbon, is_alkoxy_o) -> tuple[int, int] | None:
    """在 `carbon` 上找酯样 O 及其有效烷氧基侧 C；is_alkoxy_o 为模块特定谓词（analyzer 排除酸酐桥 O，acyl_halide 查形式电荷）。"""
    for n in carbon.GetNeighbors():
        if not is_alkoxy_o(n, carbon):
            continue
        alkoxy = _alkoxy_c_of(n, carbon)
        if alkoxy is not None:
            return n.GetIdx(), alkoxy
    return None


def _amide_n_rest(n, carbon) -> list:
    """返回酰胺 N 上除 H 与羰基碳外的其他邻居。"""
    return [x for x in n.GetNeighbors()
            if x.GetAtomicNum() != H and x.GetIdx() != carbon.GetIdx()]


def _amide_n_single(carbon, n) -> bool:
    """判断碳与 N 之间是否为单键。"""
    b = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
    return b is not None and b.GetBondType() == BondType.SINGLE


def _amide_n_substituent_ok(x) -> bool:
    """酰胺 N 上除羰基碳/H 外的一个取代基是否可接受：C，或不再连碳的羟基 O。 """
    if x.GetAtomicNum() == C:
        return True
    return x.GetAtomicNum() == O and not any(n.GetAtomicNum() == C for n in x.GetNeighbors())


def _amide_n_info(carbon) -> tuple[int, list[int]] | None:
    """返回 `carbon` 上酰胺 N 的 (n_idx, 邻居 C 索引列表)；环内 N 不作酰胺"""
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != N or n.IsInRing() or not _amide_n_single(carbon, n):
            continue
        o = _amide_n_rest(n, carbon)
        if len(o) <= 2 and all(_amide_n_substituent_ok(x) for x in o):
            return n.GetIdx(), [x.GetIdx() for x in o if x.GetAtomicNum() == C]
    return None