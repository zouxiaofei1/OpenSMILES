from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol


def preprocess(smiles: str) -> Mol | None:
    if not smiles or not str(smiles).strip():
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    return mol
