"""L0 SMILES 预处理：清洗字符串并解析为 RDKit 分子。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol


def preprocess(smiles: str) -> Mol | None:
    """清洗并解析 SMILES 为 RDKit 分子，空输入返回 None。"""
    if not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    return mol
