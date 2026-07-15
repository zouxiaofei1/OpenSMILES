"""L1 detection of free sulfonic acids R–SO2–OH / R–SO2–O⁻ (P-65.3)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol


def _dbl_o_nbs(sulfur) -> list:
    out = []
    for b in sulfur.GetBonds():
        if b.GetBondType() != BondType.DOUBLE:
            continue
        other = b.GetOtherAtom(sulfur)
        if other.GetAtomicNum() == 8:
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


def _s_core_ok(atom) -> bool:
    """Neutral tetracoordinate S with no H."""
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
        return False
    return atom.GetTotalDegree() == 4 and atom.GetFormalCharge() == 0


def _is_free_oh_o(oxygen, s_idx: int) -> bool:
    """S–OH: neutral O with H, single heavy neighbor S (not ester-O–C)."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetFormalCharge() != 0:
        return False
    if oxygen.GetTotalNumHs() < 1:
        return False
    heavies = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _is_anion_o(oxygen, s_idx: int) -> bool:
    """S–O⁻: charged O, degree 1, only neighbor S."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetFormalCharge() != -1:
        return False
    if oxygen.GetTotalDegree() != 1 or oxygen.GetTotalNumHs() != 0:
        return False
    heavies = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _free_o_of(atom) -> tuple | None:
    """Return (oxygen, is_anion) for free OH or O⁻ on S; else None."""
    s_idx = atom.GetIdx()
    for o in _sgl_nbs(atom, 8):
        if _is_free_oh_o(o, s_idx):
            return o, False
        if _is_anion_o(o, s_idx):
            return o, True
    return None


def _is_sulfonic_s(atom) -> bool:
    """S with two =O, one C, one free OH/O⁻ (not ester/amide/Cl/sulfone)."""
    if not _s_core_ok(atom) or len(_dbl_o_nbs(atom)) != 2:
        return False
    if len(_sgl_nbs(atom, 6)) != 1 or len(_sgl_nbs(atom, 8)) != 1:
        return False
    return _free_o_of(atom) is not None


def _pack_entry(atom, c, o_pair) -> dict:
    o, anion = o_pair
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(),
        "o_idx": o.GetIdx(), "anion": anion,
    }


def _entry_for(atom) -> dict | None:
    if not _is_sulfonic_s(atom):
        return None
    o_pair = _free_o_of(atom)
    if o_pair is None:
        return None
    return _pack_entry(atom, _sgl_nbs(atom, 6)[0], o_pair)


def sulfonic_acid_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry_for(a)) is not None]
