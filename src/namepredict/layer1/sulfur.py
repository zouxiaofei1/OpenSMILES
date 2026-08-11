"""L1 detection of sulfur oxy functional groups:
sulfoxides R–S(=O)–R' (P-63.3), sulfones R–SO2–R' (P-65.3.1.2),
sulfonic acids R–SO2–OH / R–SO2–O⁻ (P-65.3), sulfonates R–SO2–OR' (P-65.3.2),
sulfonamides R–SO2–NR′R″ (P-65.3), sulfonyl chlorides R–SO2–Cl (P-65.3).
"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, Cl, H, N, O, S


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


def _s_core_ok(atom) -> bool:
    """Neutral tetracoordinate S with no H."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    return atom.GetTotalDegree() == 4 and atom.GetFormalCharge() == 0


# --- Sulfoxide: R–S(=O)–R' ---

def _c_nbs(sulfur) -> list:
    return [n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == C]


def _is_sulfoxide_sulfur(atom) -> bool:
    """S with exactly one =O and two C neighbors (not sulfone/sulfide)."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 3:
        return False
    if len(_dbl_o_nbs(atom)) != 1:
        return False
    return len(_c_nbs(atom)) == 2


def _sulfoxide_entry(atom) -> dict:
    cs = _c_nbs(atom)
    return {"s_idx": atom.GetIdx(), "c1": cs[0].GetIdx(), "c2": cs[1].GetIdx()}


def sulfoxide_entries(mol: Mol) -> list[dict]:
    return [_sulfoxide_entry(a) for a in mol.GetAtoms() if _is_sulfoxide_sulfur(a)]


# --- Sulfone: R–SO2–R' ---

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


# --- Sulfonic acid: R–SO2–OH / R–SO2–O⁻ ---

def _is_free_oh_o(oxygen, s_idx: int) -> bool:
    """S–OH: neutral O with H, single heavy neighbor S (not ester-O–C)."""
    if oxygen.GetAtomicNum() != O or oxygen.GetFormalCharge() != 0:
        return False
    if oxygen.GetTotalNumHs() < 1:
        return False
    heavies = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _is_anion_o(oxygen, s_idx: int) -> bool:
    """S–O⁻: charged O, degree 1, only neighbor S."""
    if oxygen.GetAtomicNum() != O or oxygen.GetFormalCharge() != -1:
        return False
    if oxygen.GetTotalDegree() != 1 or oxygen.GetTotalNumHs() != 0:
        return False
    heavies = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _free_o_of(atom) -> tuple | None:
    """Return (oxygen, is_anion) for free OH or O⁻ on S; else None."""
    s_idx = atom.GetIdx()
    for o in _sgl_nbs(atom, O):
        if _is_free_oh_o(o, s_idx):
            return o, False
        if _is_anion_o(o, s_idx):
            return o, True
    return None


def _is_sulfonic_s(atom) -> bool:
    """S with two =O, one C, one free OH/O⁻ (not ester/amide/Cl/sulfone)."""
    if not _s_core_ok(atom) or len(_dbl_o_nbs(atom)) != 2:
        return False
    if len(_sgl_nbs(atom, C)) != 1 or len(_sgl_nbs(atom, O)) != 1:
        return False
    return _free_o_of(atom) is not None


def _pack_sulfonic_entry(atom, c, o_pair) -> dict:
    o, anion = o_pair
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(),
        "o_idx": o.GetIdx(), "anion": anion,
    }


def _sulfonic_entry_for(atom) -> dict | None:
    if not _is_sulfonic_s(atom):
        return None
    o_pair = _free_o_of(atom)
    if o_pair is None:
        return None
    return _pack_sulfonic_entry(atom, _sgl_nbs(atom, C)[0], o_pair)


def sulfonic_acid_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _sulfonic_entry_for(a)) is not None]


# --- Sulfonate ester: R–SO2–OR' ---

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


def _pack_sulfonate_entry(atom, c, o) -> dict | None:
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


def _sulfonate_entry_for(atom) -> dict | None:
    if not _is_sulfonate_s(atom):
        return None
    return _pack_sulfonate_entry(atom, _sgl_nbs(atom, C)[0], _sgl_nbs(atom, O)[0])


def sulfonate_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _sulfonate_entry_for(a)) is not None]


# --- Sulfonamide: R–SO2–NR′R″ ---

def _is_sulfonamide_s(atom) -> bool:
    """S with two =O, one C, one N (not sulfone/acid/ester/chloride)."""
    if atom.GetAtomicNum() != S or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 4 or atom.GetFormalCharge() != 0:
        return False
    if len(_dbl_o_nbs(atom)) != 2:
        return False
    return len(_sgl_nbs(atom, C)) == 1 and len(_sgl_nbs(atom, N)) == 1


def _n_rest_ok(n_atom, s_idx: int) -> bool:
    heavies = [n for n in n_atom.GetNeighbors() if n.GetAtomicNum() != H]
    if any(h.GetIdx() != s_idx and h.GetAtomicNum() not in (6,) for h in heavies):
        return False
    return sum(1 for h in heavies if h.GetAtomicNum() == C) <= 2


def _n_ok(n_atom, s_idx: int) -> bool:
    """N single-bonded only to S (+C/H); mono-N: ≤1 C outside S."""
    if n_atom.GetAtomicNum() != N or n_atom.GetIsAromatic():
        return False
    return (not None) and _n_rest_ok(n_atom, s_idx)


def _n_c_idxs(n_atom, s_idx: int) -> list[int]:
    return [x.GetIdx() for x in n_atom.GetNeighbors()
            if x.GetAtomicNum() == C and x.GetIdx() != s_idx]


def _pack_sulfonamide_entry(atom, c, n) -> dict:
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(),
        "n_idx": n.GetIdx(), "n_c_idxs": _n_c_idxs(n, atom.GetIdx()),
    }


def _sulfonamide_entry_for(atom) -> dict | None:
    if not _is_sulfonamide_s(atom):
        return None
    c, n = _sgl_nbs(atom, C)[0], _sgl_nbs(atom, N)[0]
    return _pack_sulfonamide_entry(atom, c, n) if _n_ok(n, atom.GetIdx()) else None


def sulfonamide_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _sulfonamide_entry_for(a)) is not None]


# --- Sulfonyl chloride: R–SO2–Cl ---

def _is_cl_leaf(atom) -> bool:
    """Terminal Cl single-bonded only to S (not C–Cl)."""
    if atom.GetAtomicNum() != Cl or atom.GetFormalCharge() != 0:
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetAtomicNum() == S


def _is_sulfonyl_chloride_s(atom) -> bool:
    """S with two =O, one C, one Cl (not amide/ester/OH/sulfone)."""
    if not _s_core_ok(atom) or len(_dbl_o_nbs(atom)) != 2:
        return False
    if len(_sgl_nbs(atom, C)) != 1:
        return False
    cls = _sgl_nbs(atom, Cl)
    return len(cls) == 1 and _is_cl_leaf(cls[0])


def _pack_sulfonyl_chloride_entry(atom) -> dict:
    c = _sgl_nbs(atom, C)[0]
    cl = _sgl_nbs(atom, Cl)[0]
    return {
        "s_idx": atom.GetIdx(), "c_attach": c.GetIdx(), "cl_idx": cl.GetIdx(),
    }


def _sulfonyl_chloride_entry_for(atom) -> dict | None:
    return _pack_sulfonyl_chloride_entry(atom) if _is_sulfonyl_chloride_s(atom) else None


def sulfonyl_chloride_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _sulfonyl_chloride_entry_for(a)) is not None]
