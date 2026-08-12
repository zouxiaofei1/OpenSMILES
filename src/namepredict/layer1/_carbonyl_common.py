"""Shared carbonyl-detection primitives for layer1 detectors.

Extracted verbatim from layer1/analyzer.py and layer1/acyl_halide.py, which
duplicated this whole block.  The two detectors still keep their own
`_is_ester_alkoxy_o` predicate (analyzer excludes anhydride-bridge O; acyl
halide checks formal charge) and pass it into the shared `_ester_alkoxy_of`.
"""
from __future__ import annotations

from rdkit.Chem import BondType

from namepredict.constants import C, H, N, O


def _is_single_c_oh(atom) -> bool:
    if atom.GetAtomicNum() != O or atom.GetTotalNumHs() < 1:
        return False
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == C]) == 1


def _dbl_o_on(bond, carbon) -> bool:
    if bond.GetBondType() != BondType.DOUBLE:
        return False
    return bond.GetOtherAtom(carbon).GetAtomicNum() == O


def _has_double_bonded_o(carbon) -> bool:
    return any(_dbl_o_on(b, carbon) for b in carbon.GetBonds())


def _is_carboxylate_o(atom) -> bool:
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != -1:
        return False
    return atom.GetTotalDegree() == 1 and atom.GetTotalNumHs() == 0


def _has_carboxylate_o_neighbor(carbon) -> bool:
    return any(_is_carboxylate_o(n) for n in carbon.GetNeighbors())


def _has_acid_o_neighbor(carbon) -> bool:
    return any(
        _is_single_c_oh(n) or _is_carboxylate_o(n) for n in carbon.GetNeighbors()
    )


def _is_anhydride_bridge_o(oxygen) -> bool:
    if oxygen.GetAtomicNum() != O or oxygen.GetTotalNumHs() != 0:
        return False
    cs = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() == C]
    if len(cs) != 2:
        return False
    return all(_has_double_bonded_o(c) and not _has_acid_o_neighbor(c) for c in cs)


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


def _ester_alkoxy_of(carbon, is_alkoxy_o) -> tuple[int, int] | None:
    """Find ester-like O on `carbon` whose neighbour C is a valid alkoxy side.

    `is_alkoxy_o(oxygen, carbonyl)` is the module-specific predicate: analyzer
    excludes anhydride-bridge O, acyl_halide checks formal charge.
    """
    for n in carbon.GetNeighbors():
        if not is_alkoxy_o(n, carbon):
            continue
        alkoxy = _alkoxy_c_of(n, carbon)
        if alkoxy is not None:
            return n.GetIdx(), alkoxy
    return None


def _amide_n_rest(n, carbon) -> list:
    return [x for x in n.GetNeighbors()
            if x.GetAtomicNum() != H and x.GetIdx() != carbon.GetIdx()]


def _amide_n_single(carbon, n) -> bool:
    b = carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
    return b is not None and b.GetBondType() == BondType.SINGLE


def _amide_n_info(carbon) -> tuple[int, list[int]] | None:
    """Return (n_idx, neighbour-C idx list) for the amide N on `carbon`."""
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != N or not _amide_n_single(carbon, n):
            continue
        o = _amide_n_rest(n, carbon)
        if len(o) <= 2 and all(x.GetAtomicNum() == C for x in o):
            return n.GetIdx(), [x.GetIdx() for x in o]
    return None


def _amide_n_of(carbon) -> int | None:
    info = _amide_n_info(carbon)
    return info[0] if info else None
