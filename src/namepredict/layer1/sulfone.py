"""L1 detection of open-chain dialkyl sulfones R–SO2–R' (P-65.3.1.2)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, O, S


def _dbl_o_nbs(sulfur) -> list:
    out = []
    for b in sulfur.GetBonds():
        if b.GetBondType() != BondType.DOUBLE:
            continue
        other = b.GetOtherAtom(sulfur)
        if other.GetAtomicNum() == O:
            out.append(other)
    return out


def _sgl_nbs(sulfur, z: int) -> list:
    out = []
    for b in sulfur.GetBonds():
        if b.GetBondType() != BondType.SINGLE:
            continue
        other = b.GetOtherAtom(sulfur)
        if other.GetAtomicNum() == z:
            out.append(other)
    return out


def _is_sulfone_s(atom) -> bool:
    """S with two =O and two C single bonds (not sulfoxide/sulfonamide/acid)."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 4 or atom.GetFormalCharge() != 0:
        return False
    if len(_dbl_o_nbs(atom)) != 2:
        return False
    return len(_sgl_nbs(atom, 6)) == 2


def _sulfone_entry(atom) -> dict:
    cs = _sgl_nbs(atom, 6)
    return {"s_idx": atom.GetIdx(), "c1": cs[0].GetIdx(), "c2": cs[1].GetIdx()}


def sulfone_entries(mol: Mol) -> list[dict]:
    return [_sulfone_entry(a) for a in mol.GetAtoms() if _is_sulfone_s(a)]
