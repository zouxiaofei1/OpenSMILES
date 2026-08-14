"""芳烃保留母体：苯甲酸、苯甲醛、苯乙酮、苯甲酸酯等。"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol


def _simple_alkoxy_n(mol: Mol, start: int, o_idx: int) -> int | None:
    """统计严格正构线性烷基的醇侧碳数（C1–C35，对齐 ester_alkyl 保留名表）；支链/环/重键/杂原子一律拒绝。"""
    n, cur, prev = 0, start, o_idx
    while cur is not None and n < 36:
        atom = mol.GetAtomWithIdx(cur)
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetIsAromatic():
            return None
        bond = mol.GetBondBetweenAtoms(cur, prev)
        if bond is None or bond.GetBondType() != BondType.SINGLE:
            return None
        nxt = None
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 1 or nb.GetIdx() == prev:
                continue
            if nb.GetAtomicNum() != 6:
                return None
            b = mol.GetBondBetweenAtoms(cur, nb.GetIdx())
            if b is None or b.GetBondType() != BondType.SINGLE:
                return None
            if nxt is not None:
                return None  # 支链
            nxt = nb.GetIdx()
        n, prev, cur = n + 1, cur, nxt
    return n


def _benzoate_alkoxy(mol: Mol, o_idx: int, ac: int) -> dict | None:
    """苯甲酸酯 O 侧：仅严格线性正构烷基给 alkoxy_n；其余由 L3 o_side 链路命名。"""
    n = _simple_alkoxy_n(mol, ac, o_idx)
    return {"alkoxy_n": n} if n is not None else None
