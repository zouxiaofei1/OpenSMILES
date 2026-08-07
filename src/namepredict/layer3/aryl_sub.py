"""Aryl arms: Ph / OPh / CH2Ph / OCH2Ph with recursive nested leaves (P-29.3).

Ph may carry ≤3 leaves (simple + nested Ph/OPh up to max depth).
Locants from attach=1 (lowest set).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import C, H

from namepredict.tools.aryl_depth2 import (
    _depth2_atoms_on,
)
from namepredict.tools.leaves.registry import (
    match_leaf_kind as _recurse_leaf_kind,
    nested_leaf_atoms as _nested_leaf_atoms,
)

from namepredict.constants import HALO_Z as _HALO

def _is_arom_c6(mol: Mol, atoms) -> bool:
    if len(atoms) != 6:
        return False
    return all(
        mol.GetAtomWithIdx(i).GetIsAromatic()
        and mol.GetAtomWithIdx(i).GetAtomicNum() == C
        for i in atoms
    )

def _c6_rings_at(mol: Mol, c_idx: int) -> list[set[int]]:
    return [
        set(r) for r in mol.GetRingInfo().AtomRings()
        if c_idx in r and _is_arom_c6(mol, r)
    ]

def _nb_outside(mol: Mol, i: int, ring: set[int]):
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != H and n.GetIdx() not in ring
    ]

def _is_terminal_halo(atom) -> bool:
    if atom.GetAtomicNum() not in _HALO:
        return False
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != H) == 1

def _is_terminal_me_leaf(mol: Mol, c_idx: int, ring_c: int) -> bool:
    """Methyl leaf on Ph: C only bonded to ring carbon (+H)."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != C or atom.IsInRing() or atom.GetIsAromatic():
        return False
    return all(n.GetIdx() == ring_c or n.GetAtomicNum() == H for n in atom.GetNeighbors())

def _side_leaf_kind(mol: Mol, nb, ring_i: int) -> str | None:
    if _is_terminal_halo(nb):
        return "halo"
    if nb.GetAtomicNum() == C and _is_terminal_me_leaf(mol, nb.GetIdx(), ring_i):
        return "me"
    return _recurse_leaf_kind(mol, nb, ring_i)

def _count_side_leaves(mol: Mol, ring: set[int], attach: int, parent: int) -> int | None:
    """Count terminal halo + methyl leaves on Ph (exclude parent link)."""
    n = 0
    for i in ring:
        for nb in _nb_outside(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if _side_leaf_kind(mol, nb, i) is None:
                return None
            n += 1
    return n

def _side_ok(mol: Mol, ring: set[int], attach: int, parent: int) -> bool:
    n = _count_side_leaves(mol, ring, attach, parent)
    return n is not None and n <= 3

def _phenyl_at(mol: Mol, attach: int, parent: int) -> set[int] | None:
    hits = _c6_rings_at(mol, attach)
    if len(hits) != 1 or parent in hits[0]:
        return None
    ring = hits[0]
    return ring if _side_ok(mol, ring, attach, parent) else None

def _ph_bridge_attach(mol: Mol, ph: set[int], bridge: int) -> int:
    for i in ph:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetIdx() == bridge:
                return i
    return min(ph)

def _arm_parent_of(mol: Mol, ph: set[int], attach: int) -> int:
    for n in mol.GetAtomWithIdx(attach).GetNeighbors():
        if n.GetAtomicNum() != H and n.GetIdx() not in ph:
            return n.GetIdx()
    return -1

def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]

def _is_open_ch2(atom) -> bool:
    return atom.GetAtomicNum() == C and not atom.GetIsAromatic() and not atom.IsInRing()

def _other_heavy(atom, parent: int) -> int | None:
    heavies = _heavies(atom)
    if len(heavies) != 2:
        return None
    ids = {h.GetIdx() for h in heavies}
    return None if parent not in ids else (ids - {parent}).pop()

def _ch2_ph_at(mol: Mol, ch2: int, parent: int) -> set[int] | None:
    atom = mol.GetAtomWithIdx(ch2)
    if not _is_open_ch2(atom):
        return None
    other = _other_heavy(atom, parent)
    return None if other is None else _phenyl_at(mol, other, ch2)

def _phenoxy_from_ether(mol: Mol, e: dict, parent: set[int]) -> dict | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in parent) == (c2 in parent):
        return None
    ring_c, outer = (c1, c2) if c1 in parent else (c2, c1)
    ph = _phenyl_at(mol, outer, o)
    if ph is None:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "ph": ph,
            "atoms": [o, *ph]}

def _benzyloxy_from_ether(mol: Mol, e: dict, parent: set[int]) -> dict | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in parent) == (c2 in parent):
        return None
    ring_c, outer = (c1, c2) if c1 in parent else (c2, c1)
    ph = _ch2_ph_at(mol, outer, o)
    if ph is None:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "ch2": outer,
            "ph": ph, "atoms": [o, outer, *ph]}

def _parent_link(mol: Mol, start: int, parent: set[int]) -> int | None:
    links = [
        n.GetIdx() for n in mol.GetAtomWithIdx(start).GetNeighbors()
        if n.GetIdx() in parent
    ]
    return links[0] if len(links) == 1 else None

def _phenyl_from_start(mol: Mol, start: int, parent: set[int]) -> dict | None:
    if start in parent:
        return None
    att = _parent_link(mol, start, parent)
    if att is None:
        return None
    ph = _phenyl_at(mol, start, att)
    return None if ph is None else {
        "attach": att, "outer_c": start, "ph": ph, "atoms": list(ph),
    }

def _benzyl_from_start(mol: Mol, start: int, parent: set[int]) -> dict | None:
    if start in parent:
        return None
    att = _parent_link(mol, start, parent)
    if att is None:
        return None
    ph = _ch2_ph_at(mol, start, att)
    if ph is None:
        return None
    return {"attach": att, "outer_c": start, "ch2": start, "ph": ph,
            "atoms": [start, *ph]}

def _ring_phenoxys(info: dict, parent: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _phenoxy_from_ether(mol, e, parent)
        if one is not None:
            out.append(one)
    return out

def _ring_benzyloxys(info: dict, parent: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _benzyloxy_from_ether(mol, e, parent)
        if one is not None:
            out.append(one)
    return out

def _side_c_starts(mol: Mol, parent: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for r in parent
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetAtomicNum() == C and n.GetIdx() not in parent
    ]

def _ring_phenyls(mol: Mol, parent: set[int]) -> list[dict]:
    out: list[dict] = []
    for s in _side_c_starts(mol, parent):
        one = _phenyl_from_start(mol, s, parent)
        if one is not None:
            out.append(one)
    return out

def _ring_benzyls(mol: Mol, parent: set[int]) -> list[dict]:
    out: list[dict] = []
    for s in _side_c_starts(mol, parent):
        one = _benzyl_from_start(mol, s, parent)
        if one is not None:
            out.append(one)
    return out

def _simple_leaf_atoms(mol: Mol, ph: set[int]) -> set[int]:
    out: set[int] = set()
    for i in ph:
        for nb in _nb_outside(mol, i, ph):
            if _is_terminal_halo(nb):
                out.add(nb.GetIdx())
            elif nb.GetAtomicNum() == C and _is_terminal_me_leaf(mol, nb.GetIdx(), i):
                out.add(nb.GetIdx())
    return out

def _halo_atoms_on(mol: Mol, ph: set[int], attach: int | None = None) -> set[int]:
    """Halo + methyl + depth-2 leaf atom indices on Ph (legacy name)."""
    return _simple_leaf_atoms(mol, ph) | _depth2_atoms_on(mol, ph, _nb_outside, attach)

def _arm_attach_parent(mol: Mol, p: dict) -> tuple[int, int]:
    if "ch2" in p:
        return _ph_bridge_attach(mol, p["ph"], p["ch2"]), p["ch2"]
    att = p.get("outer_c", p.get("attach"))
    return att, _arm_parent_of(mol, p["ph"], att)

def _one_arm_atoms(mol: Mol, p: dict) -> set[int]:
    att, parent = _arm_attach_parent(mol, p)
    nested = _nested_leaf_atoms(mol, p["ph"], att, parent)
    return set(p["atoms"]) | nested

def _aryl_atoms(info: dict, parent: set[int]) -> set[int]:
    mol: Mol = info["mol"]
    out: set[int] = set()
    arms = (
        _ring_phenoxys(info, parent) + _ring_benzyloxys(info, parent)
        + _ring_phenyls(mol, parent) + _ring_benzyls(mol, parent)
    )
    for p in arms:
        out |= _one_arm_atoms(mol, p)
    return out

def _aryl_sub_n(info: dict, parent: set[int]) -> int:
    mol: Mol = info["mol"]
    return (
        len(_ring_phenyls(mol, parent)) + len(_ring_phenoxys(info, parent))
        + len(_ring_benzyls(mol, parent)) + len(_ring_benzyloxys(info, parent))
    )

def _is_unfused_benzene_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True

def _phenyl_starts_set(mol: Mol, parent: set[int]) -> set[int]:
    ps = {p["outer_c"] for p in _ring_phenyls(mol, parent)}
    return ps | {p["outer_c"] for p in _ring_benzyls(mol, parent)}

def _arom_c6_ring_lists(mol: Mol) -> list[list[int]]:
    return [list(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]

def _unfused_c6_at(mol: Mol, c_idx: int) -> set[int] | None:
    """Sole unfused aromatic C6 containing c_idx, else None."""
    hits = [
        set(r) for r in _arom_c6_ring_lists(mol)
        if c_idx in r and _is_unfused_benzene_ring(mol, set(r))
    ]
    return hits[0] if len(hits) == 1 else None

def _exocyclic_fg_ring(mol: Mol, fg_c: int) -> set[int] | None:
    """Unfused C6 attached to exocyclic FG carbon (COOH/CHO/CN/...)."""
    nbs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(fg_c).GetNeighbors()
        if n.GetAtomicNum() == C and n.GetIsAromatic()
    ]
    return _unfused_c6_at(mol, nbs[0]) if len(nbs) == 1 else None

def _aryl_exclude(info: dict, parent: set[int]) -> set[int]:
    mol: Mol = info["mol"]
    return _aryl_atoms(info, parent) | _phenyl_starts_set(mol, parent)
