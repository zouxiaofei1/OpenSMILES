"""E/Z stereodescriptor prefixes (IUPAC P-91.2 / P-93.4)."""
from __future__ import annotations

from rdkit.Chem import BondStereo, Mol


def _stereo_tag(st) -> str:
    if st == BondStereo.STEREOE:
        return "(E)-"
    if st == BondStereo.STEREOZ:
        return "(Z)-"
    return ""


def _bond_stereo(mol: Mol | None, double_bond) -> str:
    if mol is None or not double_bond:
        return ""
    c1, c2 = double_bond
    bond = mol.GetBondBetweenAtoms(int(c1), int(c2))
    return _stereo_tag(bond.GetStereo()) if bond is not None else ""


def _ez_prefix(numbered: dict) -> str:
    parent = numbered.get("parent") or {}
    return _bond_stereo(parent.get("mol"), parent.get("double_bond"))


def _bond_min_loc(chain: list[int], pair) -> int | None:
    if not pair or pair[0] not in chain or pair[1] not in chain:
        return None
    return min(chain.index(pair[0]) + 1, chain.index(pair[1]) + 1)


def _ez_letter(tag: str) -> str:
    """'(E)-' → 'E'; empty → ''."""
    return tag[1] if len(tag) >= 3 and tag[0] == "(" else ""


def _ez_bond_part(mol, chain: list[int], bond) -> tuple[int, str] | None:
    loc = _bond_min_loc(chain, bond)
    letter = _ez_letter(_bond_stereo(mol, bond))
    return (loc, letter) if loc is not None and letter else None


def _ez_parts(mol, chain: list[int], bonds) -> list[tuple[int, str]]:
    """Collect (loc, letter) only for bonds with defined stereo (partial OK)."""
    parts = [_ez_bond_part(mol, chain, b) for b in bonds]
    return sorted((p for p in parts if p is not None), key=lambda x: x[0])


def _ez_multi_prefix(numbered: dict) -> str:
    """Multi-ene prefix from bonds that have stereo: (2E,6Z)- or (14Z)-."""
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), parent.get("chain") or []
    bonds = list(parent.get("double_bonds") or [])
    if mol is None or not bonds or not chain:
        return ""
    parts = _ez_parts(mol, chain, bonds)
    if not parts:
        return ""
    return f"({','.join(f'{loc}{let}' for loc, let in parts)})-"


def ez_for_parent(numbered: dict) -> str:
    """E/Z prefix for parent: multi-ene when double_bonds, else single bond."""
    parent = numbered.get("parent") or {}
    if parent.get("double_bonds"):
        return _ez_multi_prefix(numbered)
    return _ez_prefix(numbered)


# Backward-compatible alias used by unsat acid / alkenol paths.
_ez_for_alkenol = ez_for_parent
