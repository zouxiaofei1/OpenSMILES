"""L0 SMILES 预处理：清洗字符串并解析为 RDKit 分子。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.layer0.charge import normalize_acid_charge
from namepredict.layer0.tautomer import normalize_amide_tautomer


def preprocess(smiles: str) -> Mol | None:
    """清洗并解析 SMILES 为 RDKit 分子，空输入返回 None；解析后依次归一化
    酰胺烯醇互变异构为酮式酰胺、收敛酸性质子使负电荷落在最强酸位。"""
    if not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip(), sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=False, flagPossibleStereoCenters=True)
        mol = normalize_amide_tautomer(mol)
        mol = normalize_acid_charge(mol)
    except Exception:
        return None
    return mol
