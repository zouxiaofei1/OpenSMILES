"""芳烃保留母体：苯甲酸、苯甲醛、苯乙酮、苯甲酸酯等。"""
from __future__ import annotations

from rdkit.Chem import Mol

def _is_alkoxy_c(mol: Mol, cur: int, prev: int) -> bool:
    atom = mol.GetAtomWithIdx(cur)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or n.GetIdx() == prev:
            continue
        if z != 6:
            return False
    return True


def _alkoxy_next(mol: Mol, cur: int, prev: int) -> int | None:
    free = [
        n.GetIdx() for n in mol.GetAtomWithIdx(cur).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() != prev
    ]
    return free[0] if len(free) == 1 else None if not free else -1


def _simple_alkoxy_n(mol: Mol, start: int, o_idx: int) -> int | None:
    """统计线性正构烷基酯的醇侧（C1–C16）；拒绝支链/杂原子。"""
    n, cur, prev = 0, start, o_idx
    while cur is not None and n < 17:
        if not _is_alkoxy_c(mol, cur, prev) or (nxt := _alkoxy_next(mol, cur, prev)) == -1:
            return None
        n, prev, cur = n + 1, cur, nxt
    return n

#
def _benzoate_alkoxy(mol: Mol, o_idx: int, ac: int) -> dict | None:
    """线性 C1–C16 或 P-65.6 特殊烷氧基（苄基/tBu/iPr/Ph）。"""
    from namepredict.tools.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(mol, o_idx, ac)
    if side.get("alkoxy_en"):
        return side
    n = _simple_alkoxy_n(mol, ac, o_idx)
    return {**side, "alkoxy_n": n} if n is not None else None

#
def _benzoate_side_fields(side: dict) -> dict:
    return {
        "alkoxy_n": side.get("alkoxy_n"),
        "alkoxy_en": side.get("alkoxy_en") or "",
        "alkoxy_zh": side.get("alkoxy_zh") or "",
    }

