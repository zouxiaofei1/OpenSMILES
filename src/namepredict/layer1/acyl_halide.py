"""L1 detection of acyl halide R–C(=O)–X (X = Cl/Br) per IUPAC P-65.5."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import Br, C, Cl, O
from namepredict.layer1._carbonyl_common import (
    _alkoxy_c_of,
    _amide_n_of,
    _amide_n_rest,
    _amide_n_single,
    _dbl_o_on,
    _has_acid_o_neighbor,
    _has_double_bonded_o,
    _is_carboxylate_o,
    _is_single_c_oh,
)

# acyl halide detection only covers Cl/Br (P-65.5); F/I are not handled downstream
_HAL_Z = frozenset({Cl, Br})

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

def acyl_hal_of(carbon) -> tuple[int, int] | None:
    return _acyl_hal_of(carbon)
