"""L0 SMILES 预处理：清洗字符串并解析为 RDKit 分子。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.layer0.charge import normalize_acid_charge
from namepredict.layer0.tautomer import normalize_amide_tautomer


def _strip_isotopes(mol: Mol) -> Mol:
    """清除同位素标记：命名管线不产出同位素名，保留会让锚定键匹配失败。"""
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
        mol = normalize_amide_tautomer(mol)  # 电荷重定位把酰胺 O⁻ 变成中性 C(OH)=N（弱酸位受体），此时才出现可归一的酰胺烯醇位，须再跑一次（否则留下 1-hydroxyethylideneamino 式亚胺醇名，gold 取酰胺式）。
    except Exception:
        return None
    return mol
