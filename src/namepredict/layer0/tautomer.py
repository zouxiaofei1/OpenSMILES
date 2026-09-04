"""L0 酰胺烯醇互变异构归一化：非芳香中性 C(OH)=N 位点 → C(=O)-NH。

只改键级并靠 RDKit 隐氢重算完成质子迁移，不增删重原子、不改原子序；
带电 N / O⁻ 阴离子 / 显式 H / 硫类似物位点一律跳过（保守，不做质子化/阴离子改写）。
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from namepredict.constants import C, N, O


def _is_amide_enol_o(atom, carbon) -> bool:
    """判断 O 是否为可与碳上 =N 互变的烯醇羟基氧（单键、中性、仅隐氢）。"""
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != 0:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    # 只处理隐氢羟基（数据中即此形态）；显式 [OH]/H 原子跳过往 H 记账复杂化。
    return atom.GetNumImplicitHs() >= 1 and atom.GetNumExplicitHs() == 0


def _is_amide_enol_n(atom, carbon) -> bool:
    """判断 N 是否为可与碳上 -OH 互变的亚胺氮（双键、非芳香、中性、无显式 H）。"""
    if atom.GetAtomicNum() != N or atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
        return False
    if atom.GetNumExplicitHs() != 0:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    return bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE


def _amide_enol_sites(mol: Mol) -> tuple[tuple[int, int, int], ...]:
    """收集所有 (c_idx, n_idx, o_idx) 酰胺烯醇位点（非芳香中性 C-OH + C=N）。"""
    sites: list[tuple[int, int, int]] = []
    for c in mol.GetAtoms():
        if c.GetAtomicNum() != C or c.GetIsAromatic():
            continue
        o = next((n for n in c.GetNeighbors() if _is_amide_enol_o(n, c)), None)
        if o is None:
            continue
        n = next((x for x in c.GetNeighbors() if _is_amide_enol_n(x, c)), None)
        if n is not None:
            sites.append((c.GetIdx(), n.GetIdx(), o.GetIdx()))
    return tuple(sites)


def normalize_amide_tautomer(mol: Mol) -> Mol:
    """将分子内酰胺烯醇 C(OH)=N 位点归一化为酮式 C(=O)-NH；无位点或消毒失败返回原 mol。"""
    sites = _amide_enol_sites(mol)
    if not sites:
        return mol
    rw = RWMol(mol)
    for c_idx, n_idx, o_idx in sites:
        rw.RemoveBond(c_idx, n_idx)
        rw.AddBond(c_idx, n_idx, Chem.BondType.SINGLE)
        rw.RemoveBond(c_idx, o_idx)
        rw.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return mol
    Chem.AssignStereochemistry(out, force=True, cleanIt=False, flagPossibleStereoCenters=True)
    return out
