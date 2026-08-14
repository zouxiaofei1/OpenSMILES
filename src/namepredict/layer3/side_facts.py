"""Layer 3 的侧链拓扑事实：Layer 2 主链选择与 Layer 3 取代基命名共用。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C


def carbon_neighbors(mol: Mol, atom: int) -> list[int]:
    atom = mol.GetAtomWithIdx(atom)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == C]



