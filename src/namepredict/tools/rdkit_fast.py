"""把 RDKit 的 Mol.GetAtoms()/GetBonds() 换回索引循环.12.6%提升?
"""

from __future__ import annotations

_installed = False


def install() -> None:
    """安装补丁（幂等）；返回列表而非生成器，不改变元素与顺序。"""
    global _installed
    if _installed:
        return
    from rdkit import Chem

    def get_atoms(self):
        """一次索引循环取全部原子，替代逐项 Python 生成器。"""
        return [self.GetAtomWithIdx(i) for i in range(self.GetNumAtoms())]

    def get_bonds(self):
        """一次索引循环取全部键，替代逐项 Python 生成器。"""
        return [self.GetBondWithIdx(i) for i in range(self.GetNumBonds())]

    Chem.Mol.GetAtoms = get_atoms
    Chem.Mol.GetBonds = get_bonds
    _installed = True
