"""L1 detection of acyl halide R–C(=O)–X (X = Cl/Br) per IUPAC P-65.5."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import Br, C, Cl, H, N, O, HALO_Z as _HAL_Z


def _dbl_o_on(bond, carbon) -> bool:
    if bond.GetBondType() != BondType.DOUBLE:
        return False
    return bond.GetOtherAtom(carbon).GetAtomicNum() == O


def _has_double_bonded_o(carbon) -> bool:
    return any(_dbl_o_on(b, carbon) for b in carbon.GetBonds())


def _is_single_c_oh(atom) -> bool:
    if atom.GetAtomicNum() != O or atom.GetTotalNumHs() < 1:
        return False
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == C]) == 1


def _is_carboxylate_o(atom) -> bool:
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != -1:
        return False
    return atom.GetTotalDegree() == 1 and atom.GetTotalNumHs() == 0


def _has_acid_o_neighbor(carbon) -> bool:
    return any(
        _is_single_c_oh(n) or _is_carboxylate_o(n) for n in carbon.GetNeighbors()
    )


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


def _is_ester_alkoxy_o(oxygen, carbonyl) -> bool:
    if oxygen.GetAtomicNum() != O or oxygen.GetTotalNumHs() != 0:
        return False
    if oxygen.GetFormalCharge() != 0:
        return False
    return _alkoxy_c_of(oxygen, carbonyl) is not None


def _ester_alkoxy_of(carbon) -> tuple[int, int] | None:
    for n in carbon.GetNeighbors():
        if not _is_ester_alkoxy_o(n, carbon):
            continue
        alkoxy = _alkoxy_c_of(n, carbon)
        if alkoxy is not None:
            return n.GetIdx(), alkoxy
    return None


def _amide_n_rest(n, carbon) -> list:
    return [
        x for x in n.GetNeighbors()
        if x.GetAtomicNum() != H and x.GetIdx() != carbon.GetIdx()
    ]


def _amide_n_single(carbon, n) -> bool:
    b = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
    return b is not None and b.GetBondType() == BondType.SINGLE


def _amide_n_of(carbon) -> int | None:
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != N or not _amide_n_single(carbon, n):
            continue
        o = _amide_n_rest(n, carbon)
        if len(o) <= 2 and all(x.GetAtomicNum() == C for x in o):
            return n.GetIdx()
    return None


def _acyl_hal_of(carbon) -> tuple[int, int] | None:
    """Return (hal_idx, hal_z) for Cl/Br neighbor; else None."""
    for n in carbon.GetNeighbors():
        z = n.GetAtomicNum()
        if z in _HAL_Z:
            return n.GetIdx(), z
    return None


def _is_acyl_halide_carbon(atom) -> bool:
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom) or _ester_alkoxy_of(atom) is not None:
        return False
    if _amide_n_of(atom) is not None:
        return False
    return _acyl_hal_of(atom) is not None


def _entry(atom) -> dict:
    h = _acyl_hal_of(atom)
    assert h is not None
    hal_idx, hal_z = h
    return {
        "c_idx": atom.GetIdx(),
        "hal_idx": hal_idx,
        "hal_z": hal_z,
        "cl_idx": hal_idx,  # compat: L2/L3 filter still uses cl_idx
    }


def acyl_halide_entries(mol: Mol) -> list[dict]:
    return [_entry(a) for a in mol.GetAtoms() if _is_acyl_halide_carbon(a)]


def is_acyl_halide_carbon(atom) -> bool:
    return _is_acyl_halide_carbon(atom)


def acyl_hal_of(carbon) -> tuple[int, int] | None:
    return _acyl_hal_of(carbon)
