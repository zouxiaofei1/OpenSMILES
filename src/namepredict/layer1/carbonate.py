"""L1 detection of organic carbonate RO–C(=O)–OR′ (P-65.6; not ester/acid)."""
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


def _has_cooh_oh(atom) -> bool:
    return any(
        n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1 for n in atom.GetNeighbors()
    )


def _alkoxy_c_of(oxygen, carbon) -> int | None:
    cs = [
        x for x in oxygen.GetNeighbors()
        if x.GetAtomicNum() == 6 and x.GetIdx() != carbon.GetIdx()
    ]
    return cs[0].GetIdx() if len(cs) == 1 else None


def _is_ester_o(oxygen, carbon) -> bool:
    if oxygen.GetAtomicNum() != 8 or oxygen.GetTotalNumHs() != 0:
        return False
    return _bond_name(carbon, oxygen) == "SINGLE" and _alkoxy_c_of(oxygen, carbon) is not None


def _ester_o_arms(carbon) -> list[tuple[int, int]]:
    """Single-bonded O–C arms off carbonyl (ester-type alkoxy oxygens)."""
    arms: list[tuple[int, int]] = []
    for n in carbon.GetNeighbors():
        if not _is_ester_o(n, carbon):
            continue
        arms.append((n.GetIdx(), _alkoxy_c_of(n, carbon)))
    return arms


def _has_n_on(carbon) -> bool:
    return any(n.GetAtomicNum() == 7 for n in carbon.GetNeighbors())


def _is_carbonate_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6 or not _has_dbl_o(atom) or _has_cooh_oh(atom):
        return False
    if _has_n_on(atom):
        return False
    return len(_ester_o_arms(atom)) == 2


def _entry(atom) -> dict:
    arms = _ester_o_arms(atom)
    (o1, c1), (o2, c2) = arms[0], arms[1]
    return {
        "c_idx": atom.GetIdx(),
        "o1_idx": o1, "alkoxy1_c_idx": c1,
        "o2_idx": o2, "alkoxy2_c_idx": c2,
    }


def carbonate_entries(mol: Mol) -> list[dict]:
    return [_entry(a) for a in mol.GetAtoms() if _is_carbonate_carbon(a)]
