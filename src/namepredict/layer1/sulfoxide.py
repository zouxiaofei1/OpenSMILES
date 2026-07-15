"""L1 detection of open-chain dialkyl sulfoxides (P-63.3)."""
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


def _c_nbs(sulfur) -> list:
    return [n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6]


def _is_sulfoxide_sulfur(atom) -> bool:
    """S with exactly one =O and two C neighbors (not sulfone/sulfide)."""
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
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
