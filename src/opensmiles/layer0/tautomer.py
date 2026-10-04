"""L0 互变异构归一化：把次要互变体改成金标偏好的优势式（酮式/内酰胺/硫酮/胺式）。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from opensmiles.constants import C, N, O, S

ENOL_X = (O, S)  # 烯醇/烯硫醇式杂原子：环上羟基、巯基的 H 可迁往相邻环氮。
N_VALENCE = 3  # 中性氮价态上限：接收质子的氮度 + 氢数不许超过它。


def _is_amide_enol_x(atom, carbon) -> bool:
    """判断 O/S 是否为可与碳上 =N 互变的烯醇式羟基（单键、中性、带 ≥1 H）。"""
    if atom.GetAtomicNum() not in ENOL_X or atom.GetFormalCharge() != 0:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    if atom.GetDegree() != 1:  # 只与碳相连：排除以独立 H 原子形式成键的羟基
        return False
    return atom.GetTotalNumHs() >= 1  # 隐氢与显式 H 记账都接受：只认隐氢会漏新生成位点


def _is_amide_enol_n(atom, carbon) -> bool:
    """判断 N 是否为可与碳上 -XH 互变的亚胺/芳环氮（双键或芳香键、中性、可再容纳 1 个 H）。"""
    if atom.GetAtomicNum() != N or atom.GetFormalCharge() != 0:
        return False
    if atom.GetNumExplicitHs() != 0:  # 已有显式 H 的氮（吡咯型 [nH]）不接收迁移质子：再补 H 会超价。
        return False
    if atom.GetDegree() + atom.GetTotalNumHs() >= N_VALENCE:
        return False  # 连满的氮不能再接质子，否则消毒报价态超限
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    return bond is not None and bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.AROMATIC)


def _is_imine_n(atom, carbon) -> bool:
    """判断 N 是否为环外亚胺氮（中性、除碳外不连杂原子、可再容纳 1 个 H）。"""
    if atom.GetAtomicNum() != N or atom.GetFormalCharge() != 0:
        return False
    if atom.GetDegree() + atom.GetTotalNumHs() >= N_VALENCE:
        return False  # 连满的氮接不了质子：端位 =NH 与 =N-R 度+氢恰为 2
    if any(n.GetAtomicNum() != C for n in atom.GetNeighbors() if n.GetIdx() != carbon.GetIdx()):
        return False  # 肟/腙/磺酰胺的氮还连 O 或 N，不属脒，不迁移
    if atom.IsInRing():
        return False  # 亚胺氮自成环（环状脒，如咪唑并喹唑啉）时方向不定，不迁移
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    return bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE


def _is_amidine_n(atom, carbon) -> bool:
    """判断 N 是否为环内氨基氮（在环内、中性、带 H、单键或芳香键连碳）。"""
    if atom.GetAtomicNum() != N or atom.GetFormalCharge() != 0:
        return False
    if not atom.IsInRing() or atom.GetTotalNumHs() < 1:
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    return bond is not None and bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.AROMATIC)


def _is_hydroxyamino_n(atom, carbon) -> bool:
    """判断 N 是否为酰胺肟的单键羟基氨基氮（连 -OH、带 1 个 H，其 H 可迁往亚胺氮）。"""
    if atom.GetAtomicNum() != N or atom.GetFormalCharge() != 0:
        return False
    if atom.GetTotalNumHs() != 1 or atom.IsInRing():
        return False
    bond = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
        return False
    return any(n.GetAtomicNum() in ENOL_X and n.GetTotalNumHs() >= 1
               for n in atom.GetNeighbors() if n.GetIdx() != carbon.GetIdx())


def _amidoxime_candidates(mol: Mol) -> list[tuple[int, int, tuple[int, ...]]]:
    """列出 (中心碳, 亚胺氮, 羟基氨基氮) 三元组：非环酰胺肟 C(=N)-NH-OH。"""
    cands: list[tuple[int, int, tuple[int, ...]]] = []
    for c in mol.GetAtoms():
        if c.GetAtomicNum() != C:
            continue
        het = [n for n in c.GetNeighbors() if n.GetAtomicNum() != C and n.GetAtomicNum() != 1]
        if len(het) != 2 or any(n.GetAtomicNum() != N for n in het):  # 恰好两个氮、无其它杂原子
            continue
        exo = next((n for n in het if _is_imine_n(n, c)), None)
        atom = next((n for n in het if _is_hydroxyamino_n(n, c)), None)
        if exo is not None and atom is not None and exo.GetIdx() != atom.GetIdx():
            cands.append((c.GetIdx(), exo.GetIdx(), (atom.GetIdx(),)))
    return cands


def _assign(cands: list[tuple[int, int, tuple[int, ...]]]) -> tuple[tuple[int, int, int], ...]:
    """按候选氮数由少到多贪心配对：一个氮只服务一个位点，重复加氢会价态超限整体回退。"""
    used: set[int] = set()
    out: list[tuple[int, int, int]] = []
    for c_idx, y_idx, ns in sorted(cands, key=lambda s: len(s[2])):
        n_idx = next((n for n in ns if n not in used), None)
        if n_idx is not None:
            used.add(n_idx)
            out.append((c_idx, y_idx, n_idx))
    return tuple(out)

def _enol_candidates(mol: Mol) -> list[tuple[int, int, tuple[int, ...]]]:
    """列出 (中心碳, 羟基 O/S, 候选受体氮) 三元组。"""
    cands: list[tuple[int, int, tuple[int, ...]]] = []
    for c in mol.GetAtoms():
        if c.GetAtomicNum() != C:
            continue
        o = next((n for n in c.GetNeighbors() if _is_amide_enol_x(n, c)), None)
        if o is None:
            continue
        ns = tuple(x.GetIdx() for x in c.GetNeighbors() if _is_amide_enol_n(x, c))
        if ns:
            cands.append((c.GetIdx(), o.GetIdx(), ns))
    return cands


def _amidine_candidates(mol: Mol) -> list[tuple[int, int, tuple[int, ...]]]:
    """列出 (中心碳, 环外亚胺氮, 候选供体环氮) 三元组。"""
    cands: list[tuple[int, int, tuple[int, ...]]] = []
    for c in mol.GetAtoms():
        if c.GetAtomicNum() != C or not c.IsInRing():
            continue
        exo = next((x for x in c.GetNeighbors() if _is_imine_n(x, c)), None)
        if exo is None:
            continue
        ns = tuple(x.GetIdx() for x in c.GetNeighbors() if _is_amidine_n(x, c))
        if ns:
            cands.append((c.GetIdx(), exo.GetIdx(), ns))
    return cands


def _normalize_amide(rw: RWMol, sites: tuple[tuple[int, int, int], ...]) -> None:
    """就地执行烯醇式→酮式：X 上的 H 迁到 N，C=X 升为双键。"""
    for c_idx, x_idx, n_idx in sites:
        x = rw.GetAtomWithIdx(x_idx)
        if x.GetTotalNumHs() != 1:  # 羟基上多余/缺失的 H 无法靠重算隐氢弥补，本点位放弃
            continue
        x.SetNumExplicitHs(0)  # 先把 O 的 H 记账清零，否则 C=X 双键会让 O 超价
        x.SetNoImplicit(False)
        n = rw.GetAtomWithIdx(n_idx)
        n.SetNumExplicitHs(n.GetNumExplicitHs() + 1)  # 芳环 N 的隐氢不会自动补，必须显式加，否则 Kekulé 不匹配
        n.SetNoImplicit(False)
        rw.RemoveBond(c_idx, x_idx)
        rw.AddBond(c_idx, x_idx, Chem.BondType.DOUBLE)
        if rw.GetBondBetweenAtoms(c_idx, n_idx).GetBondType() != Chem.BondType.AROMATIC:
            rw.RemoveBond(c_idx, n_idx)
            rw.AddBond(c_idx, n_idx, Chem.BondType.SINGLE)


def _normalize_amidine(rw: RWMol, sites: tuple[tuple[int, int, int], ...]) -> None:
    """就地执行环外亚胺→胺式：环内 N 的 H 迁往环外 N，环内补双键。"""
    for c_idx, exo_idx, ring_idx in sites:
        exo = rw.GetAtomWithIdx(exo_idx)
        exo.SetNumExplicitHs(exo.GetNumExplicitHs() + 1)
        exo.SetNoImplicit(False)
        ring = rw.GetAtomWithIdx(ring_idx)
        ring.SetNumExplicitHs(max(0, ring.GetNumExplicitHs() - 1))
        ring.SetNoImplicit(False)
        rw.RemoveBond(c_idx, exo_idx)
        rw.AddBond(c_idx, exo_idx, Chem.BondType.SINGLE)
        if rw.GetBondBetweenAtoms(c_idx, ring_idx).GetBondType() != Chem.BondType.AROMATIC:
            rw.RemoveBond(c_idx, ring_idx)
            rw.AddBond(c_idx, ring_idx, Chem.BondType.DOUBLE)


def normalize_amide_tautomer(mol: Mol) -> Mol:
    """把烯醇/烯硫醇式 C(-XH)=N 与环外亚胺式 C(=N-H)-N< 归一化为酮式/胺式。"""
    amide = _assign(_enol_candidates(mol))
    amidine = tuple(s for s in _assign(_amidine_candidates(mol)) if s[0] not in {a[0] for a in amide})
    amidoxime = tuple(s for s in _assign(_amidoxime_candidates(mol)) if s[0] not in {a[0] for a in amide})
    if not amide and not amidine and not amidoxime:
        return mol
    rw = RWMol(mol)
    _normalize_amide(rw, amide)
    _normalize_amidine(rw, amidine + amidoxime)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return mol
    Chem.AssignStereochemistry(out, force=True, cleanIt=False, flagPossibleStereoCenters=True)
    return out
