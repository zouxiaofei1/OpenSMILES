"""芳基臂：Ph / OPh / CH2Ph / OCH2Ph 及递归嵌套叶子（P-29.3）。
Ph 至多携带 ≤3 个叶子（简单 + 嵌套 Ph/OPh 至最大深度）；位次从 attach=1 起（最小集）。
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C, H

def _is_arom_c6(mol: Mol, atoms) -> bool:
    if len(atoms) != 6:
        return False
    return all(
        mol.GetAtomWithIdx(i).GetIsAromatic()
        and mol.GetAtomWithIdx(i).GetAtomicNum() == C
        for i in atoms
    )

def _c6_rings_at(mol: Mol, c_idx: int) -> list[set[int]]:
    return [
        set(r) for r in mol.GetRingInfo().AtomRings()
        if c_idx in r and _is_arom_c6(mol, r)
    ]

def _nb_outside(mol: Mol, i: int, ring: set[int]):
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != H and n.GetIdx() not in ring
    ]



def _is_unfused_benzene_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


