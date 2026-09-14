"""L0 酰胺烯醇互变异构归一化。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from namepredict.constants import C, N, O


def _is_amide_enol_o(atom, carbon) -> bool:
    """判断 O 是否为可与碳上 =N 互变的烯醇羟基氧（单键、中性、带 ≥1 H）。"""
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != 0:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    if atom.GetDegree() != 1:  # 只与碳相连：排除以独立 H 原子（或其它重原子）形式成键的羟基，避免 H 记账复杂化。
        return False
    return atom.GetTotalNumHs() >= 1  # 隐氢与 [OH] 显式氢记账都接受：normalize_acid_charge 搬质子时写的是 explicit-H，只认隐氢会漏掉紧随其后新生成的酰胺烯醇位。


def _is_amide_enol_n(atom, carbon) -> bool:
    """判断 N 是否为可与碳上 -OH 互变的亚胺氮（双键/中性）。"""
    if atom.GetAtomicNum() != N or atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
        return False
    if atom.GetNumExplicitHs() != 0:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    return bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE


def _amide_enol_sites(mol: Mol) -> tuple[tuple[int, int, int], ...]:
    """收集所有 (c_idx, n_idx, o_idx) 酰胺烯醇位点。"""
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
    """把酰胺烯醇 C(OH)=N 位点归一化为酮式 C(=O)-NH。"""
    sites = _amide_enol_sites(mol)
    if not sites:
        return mol
    rw = RWMol(mol)
    for c_idx, n_idx, o_idx in sites:
        o = rw.GetAtomWithIdx(o_idx)
        if o.GetTotalNumHs() != 1:  # 羟基上多余/缺失的 H 无法靠重算隐氢弥补，本点位放弃
            continue
        o.SetNumExplicitHs(0)  # 羟基 O 的一个 H 搬到 N 上：先把 O 的 H 记账清零，否则 C=O 双键会让 O 价态超限、消毒失败整体回退
        o.SetNoImplicit(False)
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
