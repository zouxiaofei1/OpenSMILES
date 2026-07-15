"""L0 salt dissociation: alkali metal cation + single organic fragment.

Does not name; returns organic mol and metal metadata for L5.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

# atomic number → English metal name (IUPAC functional class salt)
_ALKALI_EN = {3: "lithium", 11: "sodium", 19: "potassium"}
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
    return a.GetAtomicNum() == 8 and a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 2


def _partition(frags: tuple[Mol, ...]) -> tuple[list[str], list[Mol]]:
    metals: list[str] = []
    organics: list[Mol] = []
    for f in frags:
        m = _alkali_en(f)
        if m is not None:
            metals.append(m)
        elif not _is_water(f):
            organics.append(f)
    return metals, organics


def _meta(metals: list[str]) -> dict | None:
    if not metals or len(set(metals)) != 1:
        return None
    en = metals[0]
    return {"metal": en, "metal_zh": _METAL_ZH[en], "n_metal": len(metals)}


def _from_frags(frags: tuple[Mol, ...]) -> tuple[Mol, dict] | None:
    metals, organics = _partition(frags)
    if len(organics) != 1:
        return None
    meta = _meta(metals)
    if meta is None:
        return None
    return organics[0], meta


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """Return (organic_mol, salt_meta). Empty meta if not a simple alkali salt."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return mol, {}
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
