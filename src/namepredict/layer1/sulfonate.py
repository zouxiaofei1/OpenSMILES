"""L1 detection of sulfonate esters R–SO2–OR′ (P-65.3.2)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, H, O, S


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


def _o_heavies(oxygen) -> list:
    return [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() != H]


def _o_links_s_c(heavies, s_idx: int) -> bool:
    has_s = any(h.GetIdx() == s_idx for h in heavies)
    has_c = any(h.GetAtomicNum() == C and h.GetIdx() != s_idx for h in heavies)
    return has_s and has_c


def _is_ester_o(oxygen, s_idx: int) -> bool:
    """Single O bonded to S and exactly one C; not OH / O- / N."""
    if oxygen.GetAtomicNum() != O or oxygen.GetFormalCharge() != 0:
        return False
    if oxygen.GetTotalNumHs() != 0 or oxygen.GetTotalDegree() != 2:
        return False
    heavies = _o_heavies(oxygen)
    return len(heavies) == 2 and _o_links_s_c(heavies, s_idx)


def _alkoxy_c_of(oxygen, s_idx: int) -> int | None:
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != s_idx:
            return n.GetIdx()
    return None


def _is_sulfonate_s(atom) -> bool:
    """S with two =O, one C, one ester-O (not amide/Cl/OH/salt)."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 4 or atom.GetFormalCharge() != 0:
        return False
    if len(_dbl_o_nbs(atom)) != 2:
        return False
    if len(_sgl_nbs(atom, C)) != 1:
        return False
    return len(_sgl_nbs(atom, O)) == 1


def _pack_entry(atom, c, o) -> dict | None:
    s_idx = atom.GetIdx()
    if not _is_ester_o(o, s_idx):
        return None
    alkoxy = _alkoxy_c_of(o, s_idx)
    if alkoxy is None:
        return None
    return {
        "s_idx": s_idx, "c_attach": c.GetIdx(),
        "o_idx": o.GetIdx(), "alkoxy_c_idx": alkoxy,
    }


def _entry_for(atom) -> dict | None:
    if not _is_sulfonate_s(atom):
        return None
    return _pack_entry(atom, _sgl_nbs(atom, C)[0], _sgl_nbs(atom, O)[0])


def sulfonate_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry_for(a)) is not None]
