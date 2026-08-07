"""Complex (recursive) leaf match handlers: nested phenyl / phenoxy / benzyl."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.aryl_depth2 import _nested_c6_at
from namepredict.tools.leaves import topo
from namepredict.tools.leaves.protocol import ArylLeafKind, Match, make_match
from namepredict.tools.leaves.simple import _FnHandler

_MAX_DEPTH = 3


def _ph_ok(mol: Mol, ph: set[int], attach: int, parent: int, depth: int) -> bool:
    from namepredict.tools.leaves.registry import all_outside_matched

    if depth > _MAX_DEPTH:
        return False
    return all_outside_matched(mol, ph, attach, parent, depth)


def match_phenyl(mol: Mol, nb, ring_i: int, depth: int) -> Match | None:
    ph = topo.nested_ph(mol, nb, ring_i)
    if ph is None:
        return None
    if not _ph_ok(mol, ph, nb.GetIdx(), ring_i, depth + 1):
        return None
    return make_match(
        ArylLeafKind.PHENYL, set(ph), ring_i,
        child_ring=ph, child_attach=nb.GetIdx(), child_parent=ring_i,
    )


def match_phenoxy(mol: Mol, nb, ring_i: int, depth: int) -> Match | None:
    outer = _phenoxy_outer(mol, nb, ring_i)
    if outer is None:
        return None
    ph = _nested_c6_at(mol, outer, {ring_i, nb.GetIdx()})
    if ph is None or not _ph_ok(mol, ph, outer, nb.GetIdx(), depth + 1):
        return None
    return make_match(
        ArylLeafKind.PHENOXY, set(ph) | {nb.GetIdx()}, ring_i,
        child_ring=ph, child_attach=outer, child_parent=nb.GetIdx(), o_idx=nb.GetIdx(),
    )


def _ether_outer_ids(o_atom, ring_i: int) -> int | None:
    if o_atom.GetAtomicNum() != 8 or o_atom.IsInRing():
        return None
    nbs = topo.heavies(o_atom)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def _phenoxy_outer(mol: Mol, o_atom, ring_i: int) -> int | None:
    outer = _ether_outer_ids(o_atom, ring_i)
    if outer is None:
        return None
    return outer if mol.GetAtomWithIdx(outer).GetIsAromatic() else None


def match_benzyl(mol: Mol, nb, ring_i: int, depth: int) -> Match | None:
    got = topo.match_benzyl_bridge(mol, nb, ring_i)
    if got is None:
        return None
    ch2, ph, outer = got
    if not _ph_ok(mol, ph, outer, ch2, depth + 1):
        return None
    return make_match(
        ArylLeafKind.BENZYL, set(ph) | {ch2}, ring_i,
        child_ring=ph, child_attach=outer, child_parent=ch2, ch2=ch2,
    )


COMPLEX_HANDLERS = [
    _FnHandler(match_phenyl, complex=True),
    _FnHandler(match_phenoxy, complex=True),
    _FnHandler(match_benzyl, complex=True),
]
