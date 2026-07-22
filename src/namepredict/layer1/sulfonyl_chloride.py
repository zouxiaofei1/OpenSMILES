"""L1 detection of sulfonyl chlorides R–SO2–Cl (P-65.3)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, Cl, H, O, S


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


def _is_cl_leaf(atom) -> bool:
    """Terminal Cl single-bonded only to S (not C–Cl)."""
    if atom.GetAtomicNum() != Cl or atom.GetFormalCharge() != 0:
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetAtomicNum() == S


def _s_core_ok(atom) -> bool:
    """Neutral tetracoordinate S with no H."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    return atom.GetTotalDegree() == 4 and atom.GetFormalCharge() == 0


def _is_sulfonyl_chloride_s(atom) -> bool:
    """S with two =O, one C, one Cl (not amide/ester/OH/sulfone)."""
    if not _s_core_ok(atom) or len(_dbl_o_nbs(atom)) != 2:
        return False
    if len(_sgl_nbs(atom, C)) != 1:
        return False
    cls = _sgl_nbs(atom, Cl)
    return len(cls) == 1 and _is_cl_leaf(cls[0])


def _pack_entry(atom) -> dict:
    c = _sgl_nbs(atom, C)[0]
    cl = _sgl_nbs(atom, Cl)[0]
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(), "cl_idx": cl.GetIdx(),
    }


def _entry_for(atom) -> dict | None:
    return _pack_entry(atom) if _is_sulfonyl_chloride_s(atom) else None


def sulfonyl_chloride_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry_for(a)) is not None]
