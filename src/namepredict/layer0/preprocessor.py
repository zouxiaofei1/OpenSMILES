"""L0 SMILES 预处理：清洗字符串并解析为 RDKit 分子。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol


def preprocess(smiles: str) -> Mol | None:
    """清洗并解析 SMILES 为 RDKit 分子，空输入返回 None。"""
    if not smiles or not str(smiles).strip():
        return None
    # 用 sanitize=False 再手动 SanitizeMol + AssignStereochemistry(cleanIt=False):
    # RDKit 自 2025 起默认解析会把隐式 H 的 [C@]/[C@@] 手性标记当作无效立体中心剥离
    # （chiral tag 丢失），导致取代基/糖苷等隐式 H 写法一进管线就丢手性。两步解析可
    # 保留这些手性中心，再加立体感知以恢复 E/Z 双键立体（cleanIt=False 避免再次剥离）。
    mol = Chem.MolFromSmiles(str(smiles).strip(), sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=False, flagPossibleStereoCenters=True)
    except Exception:
        return None
    return mol
