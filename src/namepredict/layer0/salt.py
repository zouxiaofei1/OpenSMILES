"""L0 salt dissociation: alkali metal cation or HCl + single organic fragment.

Does not name; returns organic mol and salt metadata for L5.
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


def _is_hcl_frag(mol: Mol) -> bool:
    """Neutral HCl (Cl with 1 H) or free chloride anion [Cl-]."""
    if mol.GetNumAtoms() != 1:
        return False
    a = mol.GetAtomWithIdx(0)
    if a.GetAtomicNum() != 17:
        return False
    if a.GetFormalCharge() == -1 and a.GetTotalNumHs() == 0:
        return True
    return a.GetFormalCharge() == 0 and a.GetTotalNumHs() == 1


def _bucket_frag(f: Mol, metals: list[str], organics: list[Mol]) -> int:
    """Classify one fragment; return 1 if HCl, else 0."""
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
    # Prefer alkali metal salt; exclusive of HCl co-counterion in first cut.
    if metals:
        if n_hcl:
            return None
        meta = _meta_metal(metals)
        return (organics[0], meta) if meta else None
    meta = _meta_hcl(n_hcl)
    return (organics[0], meta) if meta else None


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """Return (organic_mol, salt_meta). Empty meta if not a simple salt."""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return mol, {}
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
