"""L1 detection of simple phosphate monoesters and phosphonic acids (P-67)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol


def _is_p_oxo(atom) -> bool:
    if atom.GetAtomicNum() != 15:
        return False
    dbl = [b for b in atom.GetBonds() if b.GetBondType() == BondType.DOUBLE]
    return len(dbl) == 1 and dbl[0].GetOtherAtom(atom).GetAtomicNum() == 8


def _oh_like_o(atom, p_idx: int) -> bool:
    if atom.GetAtomicNum() != 8 or atom.GetIdx() == p_idx:
        return False
    if atom.GetFormalCharge() == -1 and atom.GetTotalDegree() == 1:
        return True
    return atom.GetTotalNumHs() >= 1 and atom.GetTotalDegree() <= 2


def _alkoxy_c(oxygen, p_idx: int) -> int | None:
    if oxygen.GetAtomicNum() != 8:
        return None
    cs = [n.GetIdx() for n in oxygen.GetNeighbors()
          if n.GetAtomicNum() == 6 and n.GetIdx() != p_idx]
    return cs[0] if len(cs) == 1 else None


def _class_nb(p, other, b, dbl_o, sgl_o, carbons) -> None:
    z = other.GetAtomicNum()
    if b.GetBondType() == BondType.DOUBLE and z == 8:
        dbl_o.append(other)
    elif b.GetBondType() == BondType.SINGLE and z == 8:
        sgl_o.append(other)
    elif z == 6:
        carbons.append(other)


def _p_neighbors(p) -> tuple[list, list, list]:
    dbl_o, sgl_o, carbons = [], [], []
    for b in p.GetBonds():
        _class_nb(p, b.GetOtherAtom(p), b, dbl_o, sgl_o, carbons)
    return dbl_o, sgl_o, carbons


def _split_sgl_o(sgl_o, p_idx: int) -> tuple[list, list]:
    pairs = [(_alkoxy_c(o, p_idx), o) for o in sgl_o]
    links = [(c, o) for c, o in pairs if c is not None]
    free = [o for c, o in pairs if c is None]
    return links, free


def _phos_links_ok(links, free, p_idx: int) -> bool:
    return len(links) == 1 and len(free) == 2 and all(
        _oh_like_o(o, p_idx) for o in free
    )


def _phosphate_entry(p) -> dict | None:
    if not _is_p_oxo(p):
        return None
    dbl_o, sgl_o, carbons = _p_neighbors(p)
    if len(dbl_o) != 1 or carbons or len(sgl_o) != 3:
        return None
    links, free = _split_sgl_o(sgl_o, p.GetIdx())
    if not _phos_links_ok(links, free, p.GetIdx()):
        return None
    return {"p_idx": p.GetIdx(), "alkoxy_c_idx": links[0][0], "o_idx": links[0][1].GetIdx()}


def _phosphonic_ok(dbl_o, sgl_o, carbons, p_idx: int) -> bool:
    if len(dbl_o) != 1 or len(carbons) != 1 or len(sgl_o) != 2:
        return False
    if any(_alkoxy_c(o, p_idx) is not None for o in sgl_o):
        return False
    return all(_oh_like_o(o, p_idx) for o in sgl_o)


def _phosphonic_entry(p) -> dict | None:
    if not _is_p_oxo(p):
        return None
    dbl_o, sgl_o, carbons = _p_neighbors(p)
    if not _phosphonic_ok(dbl_o, sgl_o, carbons, p.GetIdx()):
        return None
    return {"p_idx": p.GetIdx(), "c_idx": carbons[0].GetIdx()}


def phosphate_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _phosphate_entry(a))]


def phosphonic_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _phosphonic_entry(a))]
