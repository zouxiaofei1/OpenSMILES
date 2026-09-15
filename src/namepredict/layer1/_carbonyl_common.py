
from __future__ import annotations

from rdkit.Chem import BondType

from namepredict.constants import C, H, N, O


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


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    """返回氧上除羰基碳外的烷氧基碳索引。"""
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


def _ester_alkoxy_of(carbon, is_alkoxy_o) -> tuple[int, int] | None:
    """在 carbon 上找酯样 O 及其有效烷氧基侧 C。"""
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
    """返回 carbon 上酰胺 N 的 (n_idx, 邻居 C 索引列表)。"""
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != N or n.IsInRing() or not _amide_n_single(carbon, n):
            continue
        o = _amide_n_rest(n, carbon)
        if len(o) <= 2 and all(_amide_n_substituent_ok(x) for x in o):
            return n.GetIdx(), [x.GetIdx() for x in o if x.GetAtomicNum() == C]
    return None