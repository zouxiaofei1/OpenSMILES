"""L1 detection of carbamate R2N-C(=O)-OR (P-65; not ester/amide)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _bond_name(a, b) -> str:
    bond = a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
    return bond.GetBondType().name if bond is not None else ""


def _has_dbl_o(carbon) -> bool:
    return any(
        n.GetAtomicNum() == 8 and _bond_name(carbon, n) == "DOUBLE"
        for n in carbon.GetNeighbors()
    )


def _sgl_o_nbs(carbon):
    return [
        n for n in carbon.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 0
        and _bond_name(carbon, n) == "SINGLE"
    ]


def _alkoxy_of(carbon) -> tuple[int, int] | None:
    for n in _sgl_o_nbs(carbon):
        cs = [x for x in n.GetNeighbors()
              if x.GetAtomicNum() == 6 and x.GetIdx() != carbon.GetIdx()]
        if len(cs) == 1:
            return n.GetIdx(), cs[0].GetIdx()
    return None


def _n_rest(n, carbon):
    return [x for x in n.GetNeighbors()
            if x.GetAtomicNum() != 1 and x.GetIdx() != carbon.GetIdx()]


def _n_info(carbon) -> tuple[int, list[int]] | None:
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != 7 or _bond_name(carbon, n) != "SINGLE":
            continue
        rest = _n_rest(n, carbon)
        if len(rest) <= 2 and all(x.GetAtomicNum() == 6 for x in rest):
            return n.GetIdx(), [x.GetIdx() for x in rest]
    return None


def _has_cooh_oh(atom) -> bool:
    return any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1
               for n in atom.GetNeighbors())


def _is_carbamate_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_dbl_o(atom) or _has_cooh_oh(atom):
        return False
    return _alkoxy_of(atom) is not None and _n_info(atom) is not None


def _entry(atom) -> dict:
    o_idx, alkoxy_c = _alkoxy_of(atom)
    n_idx, n_cs = _n_info(atom)
    return {
        "c_idx": atom.GetIdx(), "o_idx": o_idx, "alkoxy_c_idx": alkoxy_c,
        "n_idx": n_idx, "n_c_idxs": n_cs,
    }


def carbamate_entries(mol: Mol) -> list[dict]:
    return [_entry(a) for a in mol.GetAtoms() if _is_carbamate_carbon(a)]
