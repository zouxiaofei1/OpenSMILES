"""Alkenedioic acid name assembly and E/Z stereo prefix."""
from __future__ import annotations

from rdkit.Chem import BondStereo, Mol

from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH


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


def _alkenedioic_names(n: int, locant: int | None, ez: str) -> tuple[str, str] | None:
    en = ALKANE_EN.get(n)
    zh = ALKANE_ZH.get(n)
    if not en or not zh or locant is None or n < 3:
        return None
    stem_en = en[:-3]  # butane -> but
    return (
        f"{ez}{stem_en}-{locant}-enedioic acid",
        f"{ez}{zh[0]}-{locant}-烯二酸",
    )


def alkenedioic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    return _alkenedioic_names(n, numbered.get("ene_locant"), _ez_prefix(numbered))
