"""L0 SMILES 预处理：清洗字符串并解析为 RDKit 分子。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.constants import ISO_H_PROPS, ISO_NUCLIDE_PROP
from namepredict.layer0.charge import normalize_acid_charge
from namepredict.layer0.tautomer import normalize_amide_tautomer

_H_ISO_PROP = {mass: prop for mass, prop, _, _ in ISO_H_PROPS}  # 氢同位素质量数 → 重原子属性名


def _capture_isotopes(mol: Mol) -> None:
    """把同位素记到重原子属性上（供 L5 渲染同位素名），清零前须先采集。"""
    for atom in mol.GetAtoms():
        iso = atom.GetIsotope()
        if not iso:
            continue
        if atom.GetAtomicNum() != 1:  # 非氢核素：记在核素原子自身，如 "13C"/"18F"
            atom.SetProp(ISO_NUCLIDE_PROP, f"{iso}{atom.GetSymbol()}")
            continue
        key = _H_ISO_PROP.get(iso)  # 氘/氚：计数记到其唯一重原子邻居上
        nb = atom.GetNeighbors()
        if key is None or len(nb) != 1:
            continue
        heavy = nb[0]
        heavy.SetIntProp(key, (heavy.GetIntProp(key) if heavy.HasProp(key) else 0) + 1)


def _strip_isotopes(mol: Mol) -> Mol:
    """清除同位素标记：命名管线用重原子属性渲染同位素名，保留会让锚定键匹配失败。"""
    _capture_isotopes(mol)
    for atom in mol.GetAtoms():
        if atom.GetIsotope():
            atom.SetIsotope(0)
    return mol


def preprocess(smiles: str) -> Mol | None:
    """清洗解析 SMILES 为 RDKit 分子；空输入/失败返回 None。"""
    if not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip(), sanitize=False)
    if mol is None:
        return None
    mol = _strip_isotopes(mol)
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=False, flagPossibleStereoCenters=True)
        mol = normalize_amide_tautomer(mol)
        mol = normalize_acid_charge(mol)
        mol = normalize_amide_tautomer(mol)  # 电荷重定位后才出现可归一的酰胺烯醇位，须再跑一次
    except Exception:
        return None
    return mol
