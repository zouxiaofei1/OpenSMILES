"""L1 detection of sulfonamides R–SO2–NR′R″ (P-65.3)."""
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


def _is_sulfonamide_s(atom) -> bool:
    """S with two =O, one C, one N (not sulfone/acid/ester/chloride)."""
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 4 or atom.GetFormalCharge() != 0:
        return False
    if len(_dbl_o_nbs(atom)) != 2:
        return False
    return len(_sgl_nbs(atom, 6)) == 1 and len(_sgl_nbs(atom, 7)) == 1


def _n_ok(n_atom, s_idx: int) -> bool:
    """N single-bonded only to S (+C/H); mono-N: ≤1 C outside S."""
    if n_atom.GetAtomicNum() != 7 or n_atom.GetIsAromatic():
        return False
    heavies = [n for n in n_atom.GetNeighbors() if n.GetAtomicNum() != 1]
    if any(h.GetIdx() != s_idx and h.GetAtomicNum() not in (6,) for h in heavies):
        return False
    n_c = sum(1 for h in heavies if h.GetAtomicNum() == 6)
    return n_c <= 1


def _n_c_idxs(n_atom, s_idx: int) -> list[int]:
    return [x.GetIdx() for x in n_atom.GetNeighbors()
            if x.GetAtomicNum() == 6 and x.GetIdx() != s_idx]


def _pack_entry(atom, c, n) -> dict:
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(),
        "n_idx": n.GetIdx(), "n_c_idxs": _n_c_idxs(n, atom.GetIdx()),
    }


def _entry_for(atom) -> dict | None:
    if not _is_sulfonamide_s(atom):
        return None
    c, n = _sgl_nbs(atom, 6)[0], _sgl_nbs(atom, 7)[0]
    return _pack_entry(atom, c, n) if _n_ok(n, atom.GetIdx()) else None


def sulfonamide_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry_for(a)) is not None]
