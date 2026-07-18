"""Terminal parent ownership: immutable owned_atoms finalization."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _dbl_o_idx


def _chain_atoms(parent: dict) -> set[int]:
    return set(parent.get("chain") or [])


def _add_opt(out: set[int], idx: int | None) -> set[int]:
    if idx is not None:
        out.add(idx)
    return out


def _amide_n_from_c(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _amide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    c_idx = parent.get("amide_c_idx")
    if c_idx is None:
        return set()
    out = {c_idx}
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _amide_n_from_c(mol, c_idx))


def _single_o_idx(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != "SINGLE":
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == 8:
            return other.GetIdx()
    return None


def _acid_o_atoms(mol: Mol, c_idx: int) -> set[int]:
    out: set[int] = set()
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _single_o_idx(mol, c_idx))


def _cooh_c_idxs(parent: dict) -> list[int]:
    multi = parent.get("cooh_c_idxs")
    if multi:
        return [int(x) for x in multi]
    one = parent.get("cooh_c_idx")
    return [int(one)] if one is not None else []


def _acid_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    out: set[int] = set()
    for c in _cooh_c_idxs(parent):
        out.add(c)
        out |= _acid_o_atoms(mol, c)
    return out


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    kind = parent.get("kind")
    out: set[int] = set()
    if kind in ("benzamide", "amide"):
        out |= _amide_fg_atoms(mol, parent)
    out |= _acid_fg_atoms(mol, parent)
    return out


def compute_owned_atoms(parent: dict, mol: Mol) -> frozenset[int]:
    """Union chain + kind-specific FG atoms (terminal ownership set)."""
    return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))


def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
    """Copy candidate once with immutable owned_atoms frozenset."""
    if isinstance(parent.get("owned_atoms"), frozenset):
        return parent
    return {**parent, "owned_atoms": compute_owned_atoms(parent, mol)}
