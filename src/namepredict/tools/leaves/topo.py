"""Shared topology helpers for leaf handlers."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.aryl_depth2 import (
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
from namepredict.tools.leaves.protocol import ArylLeafKind, make_match

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
    return make_match(ArylLeafKind.HALOGEN, {nb.GetIdx()}, ring_i, z=nb.GetAtomicNum())


def match_me(mol: Mol, nb, ring_i: int, depth: int):
    if nb.GetAtomicNum() != 6 or not is_me_leaf(mol, nb.GetIdx(), ring_i):
        return None
    return make_match(ArylLeafKind.METHYL, {nb.GetIdx()}, ring_i)


def match_cf3(mol: Mol, nb, ring_i: int, depth: int):
    if nb.GetAtomicNum() != 6 or not _is_cf3_leaf(mol, nb.GetIdx(), ring_i):
        return None
    return make_match(ArylLeafKind.TRIFLUOROMETHYL, _leaf_atoms_cf3(mol, nb.GetIdx()), ring_i)


def match_nitro(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_nitro(nb):
        return None
    return make_match(ArylLeafKind.NITRO, _leaf_atoms_nitro(nb), ring_i)


def match_hydroxy(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_oh(nb, ring_i):
        return None
    return make_match(ArylLeafKind.HYDROXY, {nb.GetIdx()}, ring_i)


def match_amino(mol: Mol, nb, ring_i: int, depth: int):
    if not _is_terminal_nh2(nb, ring_i):
        return None
    return make_match(ArylLeafKind.AMINO, {nb.GetIdx()}, ring_i)


def match_alkoxy(mol: Mol, nb, ring_i: int, depth: int):
    n = _alkoxy_n(mol, nb, ring_i)
    if n is None or n not in _ALKOXY_EN:
        return None
    atoms = _leaf_atoms_alkoxy(mol, nb.GetIdx(), ring_i)
    return make_match(ArylLeafKind.ALKOXY, atoms, ring_i, n=n)


def _open_c_only(mol: Mol, idx: int, prev: int) -> bool:
    a = mol.GetAtomWithIdx(idx)
    if a.GetAtomicNum() != 6 or a.IsInRing() or a.GetIsAromatic():
        return False
    return all(
        n.GetAtomicNum() in (1, 6) and (n.GetIdx() == prev or not n.GetIsAromatic())
        for n in a.GetNeighbors()
    )


def _one_fwd_c(mol: Mol, cur: int, prev: int) -> int | None:
    nxts = [
        n.GetIdx() for n in mol.GetAtomWithIdx(cur).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() != prev
    ]
    return nxts[0] if len(nxts) == 1 else None


def _n_alkyl_step(mol: Mol, path: list[int], prev: int) -> int | None:
    nxt = _one_fwd_c(mol, path[-1], prev)
    if nxt is None or not _open_c_only(mol, nxt, path[-1]):
        return None
    return nxt


def _n_alkyl_path(mol: Mol, start: int, ring_i: int) -> list[int] | None:
    """Linear open C2–C4 chain from ring carbon (exclude methyl C1)."""
    if not _open_c_only(mol, start, ring_i):
        return None
    path, prev = [start], ring_i
    while len(path) < 4 and (nxt := _n_alkyl_step(mol, path, prev)) is not None:
        prev = path[-1]
        path.append(nxt)
    return path if 2 <= len(path) <= 4 else None


def match_n_alkyl(mol: Mol, nb, ring_i: int, depth: int):
    if nb.GetAtomicNum() != 6:
        return None
    path = _n_alkyl_path(mol, nb.GetIdx(), ring_i)
    if path is None:
        return None
    return make_match(ArylLeafKind.N_ALKYL, set(path), ring_i, n=len(path))


def _s_me_outer(mol: Mol, s_atom, ring_i: int) -> int | None:
    if s_atom.GetAtomicNum() != 16 or s_atom.IsInRing():
        return None
    nbs = heavies(s_atom)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def match_methylthio(mol: Mol, nb, ring_i: int, depth: int):
    """Ring–S–Me terminal."""
    me = _s_me_outer(mol, nb, ring_i)
    if me is None or not is_me_leaf(mol, me, nb.GetIdx()):
        return None
    return make_match(ArylLeafKind.METHYLSULFANYL, {nb.GetIdx(), me}, ring_i)


def _cyano_n(nb, ring_i: int):
    if nb.GetAtomicNum() != 6 or nb.IsInRing() or nb.GetIsAromatic():
        return None
    nbs = heavies(nb)
    if len(nbs) != 2 or {n.GetAtomicNum() for n in nbs} != {6, 7}:
        return None
    n_atom = next(n for n in nbs if n.GetAtomicNum() == 7)
    c_link = next(n for n in nbs if n.GetAtomicNum() == 6)
    if c_link.GetIdx() != ring_i or len(heavies(n_atom)) != 1:
        return None
    return n_atom


def match_cyano(mol: Mol, nb, ring_i: int, depth: int):
    """Ring–C≡N (nitrile carbon attached to ring)."""
    n_atom = _cyano_n(nb, ring_i)
    if n_atom is None:
        return None
    return make_match(ArylLeafKind.CYANO, {nb.GetIdx(), n_atom.GetIdx()}, ring_i)


def nested_ph(mol: Mol, nb, ring_i: int) -> set[int] | None:
    if nb.GetAtomicNum() != 6 or not nb.GetIsAromatic():
        return None
    ph = _nested_c6_at(mol, nb.GetIdx(), {ring_i})
    return None if ph is None or ring_i in ph else ph


def _ch2_other(nb, ring_i: int) -> int | None:
    if nb.GetAtomicNum() != 6 or nb.IsInRing() or nb.GetIsAromatic():
        return None
    nbs = heavies(nb)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def match_benzyl_bridge(mol: Mol, nb, ring_i: int) -> tuple[int, set[int], int] | None:
    """Ring–CH2–Ph: return (ch2_idx, ph_set, ph_attach) if topology matches."""
    other = _ch2_other(nb, ring_i)
    if other is None:
        return None
    ph = _nested_c6_at(mol, other, {ring_i, nb.GetIdx()})
    if ph is None or other not in ph:
        return None
    return nb.GetIdx(), ph, other
