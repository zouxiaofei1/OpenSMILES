"""L0 盐解离：碱金属或 HCl 盐；返回有机 mol 与盐元数据供 L5。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.constants import ALKALI_EN, Cl, METAL_ZH


def _alkali_en(mol: Mol) -> str | None:
    """若 mol 为单一 +1 碱金属原子则返回其英文名。"""
    if mol.GetNumAtoms() != 1:
        return None
    atom = mol.GetAtomWithIdx(0)
    if atom.GetFormalCharge() != 1:
        return None
    return ALKALI_EN.get(atom.GetAtomicNum())

def _is_hcl_frag(mol: Mol) -> bool:
    """中性 HCl（带 1 个 H 的 Cl）或游离氯阴离子 [Cl-]。"""
    if mol.GetNumAtoms() != 1:
        return False
    a = mol.GetAtomWithIdx(0)
    if a.GetAtomicNum() != Cl:
        return False
    if a.GetFormalCharge() == -1 and a.GetTotalNumHs() == 0:
        return True
    return a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 1


def _from_frags(frags: tuple[Mol, ...]) -> tuple[Mol, dict] | None:
    """从片段中提取单一有机分子与盐元数据；不满足则返回 None。"""
    metals: list[str] = []
    organics: list[Mol] = []
    n_hcl = 0
    for f in frags:
        m = _alkali_en(f)
        if m is not None:
            metals.append(m)
        elif _is_hcl_frag(f):
            n_hcl += 1
        else:
            organics.append(f)
    if len(organics) != 1:
        return None
    if metals:  # 优先碱金属盐；第一版排除 HCl 共抗衡离子，且金属种类须唯一。
        if n_hcl or len(set(metals)) != 1:
            return None
        meta = {"metal": metals[0], "metal_zh": METAL_ZH[metals[0]], "n_metal": len(metals)}
    else:
        if n_hcl != 1:
            return None
        meta = {"acid_salt": "hydrochloride", "acid_salt_zh": "盐酸盐"}
    return organics[0], meta


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """返回有机 mol 与盐元数据；非简单盐时 meta 为空。"""
    if len(Chem.GetMolFrags(mol)) < 2:
        return mol, {}
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
