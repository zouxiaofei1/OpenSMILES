"""L0 盐解离：碱金属阳离子或 HCl + 单一有机片段；不命名，返回有机 mol 与盐元数据供 L5 使用。
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.constants import Cl, K, Li, Na, O

# 原子序数 → 英文金属名（IUPAC 官能团类盐）
_ALKALI_EN = {Li: "lithium", Na: "sodium", K: "potassium"}
_METAL_ZH = {"lithium": "锂", "sodium": "钠", "potassium": "钾"}


def _alkali_en(mol: Mol) -> str | None:
    if mol.GetNumAtoms() != 1:
        return None
    atom = mol.GetAtomWithIdx(0)
    if atom.GetFormalCharge() != 1:
        return None
    return _ALKALI_EN.get(atom.GetAtomicNum())


def _is_water(mol: Mol) -> bool:
    if mol.GetNumAtoms() != 1:
        return False
    a = mol.GetAtomWithIdx(0)
    return a.GetAtomicNum() == O and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 2


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


def _bucket_frag(f: Mol, metals: list[str], organics: list[Mol]) -> int:
    """对一个片段分类；若为 HCl 返回 1，否则返回 0。"""
    m = _alkali_en(f)
    if m is not None:
        metals.append(m)
        return 0
    if _is_hcl_frag(f):
        return 1
    if not _is_water(f):
        organics.append(f)
    return 0


def _partition(frags: tuple[Mol, ...]) -> tuple[list[str], list[Mol], int]:
    metals: list[str] = []
    organics: list[Mol] = []
    n_hcl = sum(_bucket_frag(f, metals, organics) for f in frags)
    return metals, organics, n_hcl


def _meta_metal(metals: list[str]) -> dict | None:
    if not metals or len(set(metals)) != 1:
        return None
    en = metals[0]
    return {"metal": en, "metal_zh": _METAL_ZH[en], "n_metal": len(metals)}


def _meta_hcl(n_hcl: int) -> dict | None:
    if n_hcl != 1:
        return None
    return {"acid_salt": "hydrochloride", "acid_salt_zh": "盐酸盐"}


def _from_frags(frags: tuple[Mol, ...]) -> tuple[Mol, dict] | None:
    metals, organics, n_hcl = _partition(frags)
    if len(organics) != 1:
        return None
    # 优先碱金属盐；第一版排除 HCl 共抗衡离子。
    if metals:
        if n_hcl:
            return None
        meta = _meta_metal(metals)
        return (organics[0], meta) if meta else None
    meta = _meta_hcl(n_hcl)
    return (organics[0], meta) if meta else None


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """返回 (organic_mol, salt_meta)；非简单盐时 meta 为空。"""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return mol, {}
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
