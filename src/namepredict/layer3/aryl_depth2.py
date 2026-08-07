"""Depth-2 simple leaves on Ph arms: alkoxy / nitro / CF3 (P-29.3).

Topology + prefix only; ring walk / attach locants stay in aryl_sub.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.side_alkoxy import _linear_alkyl_atoms
from namepredict.layer3.side_alkyl import _is_cf3_carbon

_ALKOXY_EN = {1: "methoxy", 2: "ethoxy", 3: "propoxy", 4: "butoxy"}
_ALKOXY_ZH = {1: "甲氧基", 2: "乙氧基", 3: "丙氧基", 4: "丁氧基"}


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _nitro_os_ok(nbs: list) -> bool:
    os = [n for n in nbs if n.GetAtomicNum() == 8]
    cs = [n for n in nbs if n.GetAtomicNum() == 6]
    if len(os) != 2 or len(cs) != 1:
        return False
    return all(len(_heavies(o)) == 1 for o in os)


def _is_terminal_nitro(atom) -> bool:
    """Ring-NO2: N with exactly two O neighbors (no further heavy)."""
    if atom.GetAtomicNum() != 7:
        return False
    nbs = _heavies(atom)
    return len(nbs) == 3 and _nitro_os_ok(nbs)


def _ether_outer(o_atom, ring_i: int) -> int | None:
    if o_atom.GetAtomicNum() != 8 or o_atom.IsInRing():
        return None
    nbs = _heavies(o_atom)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def _alkoxy_n(mol: Mol, o_atom, ring_i: int) -> int | None:
    """Ether O on ring → linear C1–C4 alkyl length, else None."""
    outer = _ether_outer(o_atom, ring_i)
    if outer is None:
        return None
    path = _linear_alkyl_atoms(mol, outer, o_atom.GetIdx())
    return len(path) if path is not None and 1 <= len(path) <= 4 else None


def _is_cf3_leaf(mol: Mol, c_idx: int, ring_i: int) -> bool:
    if not _is_cf3_carbon(mol, c_idx):
        return False
    atom = mol.GetAtomWithIdx(c_idx)
    cs = [n for n in _heavies(atom) if n.GetAtomicNum() == 6]
    return len(cs) == 1 and cs[0].GetIdx() == ring_i


def _alkoxy_kind(n: int | None) -> str | None:
    return _ALKOXY_EN.get(n) if n else None


def _is_terminal_oh(nb, ring_i: int) -> bool:
    """Phenolic OH: O bonded only to ring carbon (+H)."""
    if nb.GetAtomicNum() != 8 or nb.IsInRing():
        return False
    nbs = _heavies(nb)
    return len(nbs) == 1 and nbs[0].GetIdx() == ring_i


def _is_terminal_nh2(nb, ring_i: int) -> bool:
    """Primary amino on ring: N bonded only to ring carbon (+H)."""
    if nb.GetAtomicNum() != 7 or nb.IsInRing():
        return False
    nbs = _heavies(nb)
    return len(nbs) == 1 and nbs[0].GetIdx() == ring_i


def _is_arom_c6_set(mol: Mol, atoms: set[int]) -> bool:
    if len(atoms) != 6:
        return False
    return all(
        mol.GetAtomWithIdx(i).GetIsAromatic()
        and mol.GetAtomWithIdx(i).GetAtomicNum() == 6
        for i in atoms
    )


def _nested_c6_at(mol: Mol, c_idx: int, parent_ring: set[int]) -> set[int] | None:
    """Sole aromatic C6 at c_idx that does not contain parent_ring atoms."""
    hits = [
        set(r) for r in mol.GetRingInfo().AtomRings()
        if c_idx in r and _is_arom_c6_set(mol, set(r)) and not (set(r) & parent_ring)
    ]
    return hits[0] if len(hits) == 1 else None


def _nested_ph_unsub(mol: Mol, ph: set[int], attach: int, parent_ring_c: int) -> bool:
    """True if nested Ph has no outside heavy except parent link."""
    for i in ph:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            if nb.GetAtomicNum() == 1 or nb.GetIdx() in ph:
                continue
            if i == attach and nb.GetIdx() == parent_ring_c:
                continue
            return False
    return True


def _is_nested_phenyl(mol: Mol, nb, ring_i: int) -> bool:
    """Unsubstituted Ph leaf attached to arm ring carbon ring_i."""
    if nb.GetAtomicNum() != 6 or not nb.GetIsAromatic():
        return False
    ph = _nested_c6_at(mol, nb.GetIdx(), {ring_i})
    if ph is None or ring_i in ph:
        return False
    return _nested_ph_unsub(mol, ph, nb.GetIdx(), ring_i)


def _depth2_hetero(mol: Mol, nb, ring_i: int) -> str | None:
    if _is_terminal_nitro(nb):
        return "nitro"
    if _is_terminal_oh(nb, ring_i):
        return "hydroxy"
    if _is_terminal_nh2(nb, ring_i):
        return "amino"
    return _alkoxy_kind(_alkoxy_n(mol, nb, ring_i))


def _depth2_kind(mol: Mol, nb, ring_i: int) -> str | None:
    """nitro | cf3 | alkoxy | hydroxy | amino | phenyl leaf kind."""
    if nb.GetAtomicNum() == 6 and _is_cf3_leaf(mol, nb.GetIdx(), ring_i):
        return "cf3"
    if _is_nested_phenyl(mol, nb, ring_i):
        return "phenyl"
    return _depth2_hetero(mol, nb, ring_i)


def _leaf_atoms_nitro(nb) -> set[int]:
    return {nb.GetIdx()} | {o.GetIdx() for o in _heavies(nb) if o.GetAtomicNum() == 8}


def _leaf_atoms_cf3(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {c_idx} | {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 9}


def _leaf_atoms_alkoxy(mol: Mol, o_idx: int, ring_i: int) -> set[int]:
    o = mol.GetAtomWithIdx(o_idx)
    nbs = _heavies(o)
    outer = next(n.GetIdx() for n in nbs if n.GetIdx() != ring_i)
    path = _linear_alkyl_atoms(mol, outer, o_idx) or []
    return {o_idx, *path}


def _leaf_atoms_nested_ph(mol: Mol, nb, ring_i: int) -> set[int]:
    ph = _nested_c6_at(mol, nb.GetIdx(), {ring_i})
    return set(ph) if ph is not None else {nb.GetIdx()}


def _leaf_atoms_one(mol: Mol, nb, ring_i: int, skip_attach: int | None = None) -> set[int]:
    if skip_attach is not None and ring_i == skip_attach and _is_nested_phenyl(mol, nb, ring_i):
        return set()
    kind = _depth2_kind(mol, nb, ring_i)
    if kind == "nitro":
        return _leaf_atoms_nitro(nb)
    if kind == "cf3":
        return _leaf_atoms_cf3(mol, nb.GetIdx())
    if kind == "phenyl":
        return _leaf_atoms_nested_ph(mol, nb, ring_i)
    return _leaf_atoms_simple(mol, nb, ring_i, kind)


def _leaf_atoms_simple(mol: Mol, nb, ring_i: int, kind: str | None) -> set[int]:
    if kind in _ALKOXY_EN.values():
        return _leaf_atoms_alkoxy(mol, nb.GetIdx(), ring_i)
    if kind in ("hydroxy", "amino"):
        return {nb.GetIdx()}
    return set()


def _depth2_atoms_on(mol: Mol, ph: set[int], nb_outside, attach: int | None = None) -> set[int]:
    out: set[int] = set()
    for i in ph:
        for nb in nb_outside(mol, i, ph):
            out |= _leaf_atoms_one(mol, nb, i, attach)
    return out
