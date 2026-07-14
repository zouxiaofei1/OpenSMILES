"""Shared topology helpers for leaf handlers."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import (
    _ALKOXY_EN,
    _alkoxy_n,
    _heavies,
    _is_cf3_leaf,
    _is_terminal_nh2,
    _is_terminal_nitro,
    _is_terminal_oh,
    _leaf_atoms_alkoxy,
    _leaf_atoms_cf3,
    _leaf_atoms_nitro,
    _nested_c6_at,
)
from namepredict.layer2.leaves.protocol import make_match

_HALO = frozenset({9, 17, 35, 53})


def heavies(atom) -> list:
    return _heavies(atom)


def is_halo(nb) -> bool:
    if nb.GetAtomicNum() not in _HALO:
        return False
    return sum(1 for n in nb.GetNeighbors() if n.GetAtomicNum() != 1) == 1


def is_me_leaf(mol: Mol, c_idx: int, ring_c: int) -> bool:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetIsAromatic():
        return False
    return all(n.GetIdx() == ring_c or n.GetAtomicNum() == 1 for n in atom.GetNeighbors())


def nb_out(mol: Mol, i: int, ring: set[int]) -> list:
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ring
    ]


def match_halo(mol: Mol, nb, ring_i: int, depth: int):
    if not is_halo(nb):
        return None
    return make_match("halo", {nb.GetIdx()}, ring_i, z=nb.GetAtomicNum())


def match_me(mol: Mol, nb, ring_i: int, depth: int):
    if nb.GetAtomicNum() != 6 or not is_me_leaf(mol, nb.GetIdx(), ring_i):
        return None
    return make_match("me", {nb.GetIdx()}, ring_i)


def match_cf3(mol: Mol, nb, ring_i: int, depth: int):
    if nb.GetAtomicNum() != 6 or not _is_cf3_leaf(mol, nb.GetIdx(), ring_i):
        return None
    return make_match("cf3", _leaf_atoms_cf3(mol, nb.GetIdx()), ring_i)


def match_nitro(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_nitro(nb):
        return None
    return make_match("nitro", _leaf_atoms_nitro(nb), ring_i)


def match_hydroxy(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_oh(nb, ring_i):
        return None
    return make_match("hydroxy", {nb.GetIdx()}, ring_i)


def match_amino(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_nh2(nb, ring_i):
        return None
    return make_match("amino", {nb.GetIdx()}, ring_i)


def match_alkoxy(mol: Mol, nb, ring_i: int, depth: int):
    n = _alkoxy_n(mol, nb, ring_i)
    if n is None or n not in _ALKOXY_EN:
        return None
    atoms = _leaf_atoms_alkoxy(mol, nb.GetIdx(), ring_i)
    return make_match("alkoxy", atoms, ring_i, n=n)


def nested_ph(mol: Mol, nb, ring_i: int) -> set[int] | None:
    if nb.GetAtomicNum() != 6 or not nb.GetIsAromatic():
        return None
    ph = _nested_c6_at(mol, nb.GetIdx(), {ring_i})
    return None if ph is None or ring_i in ph else ph
