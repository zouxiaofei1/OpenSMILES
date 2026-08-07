"""Leaf registry: complex handlers first, then simple."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.tools.leaves.complex_h import COMPLEX_HANDLERS
from namepredict.tools.leaves.protocol import ArylLeafKind, LeafTopology, Match
from namepredict.tools.leaves.simple import SIMPLE_HANDLERS
from namepredict.tools.leaves.topo import nb_out

# complex before simple (OPh must win over future ether leaves)
_HANDLERS = list(COMPLEX_HANDLERS) + list(SIMPLE_HANDLERS)


def match_leaf(mol: Mol, nb, ring_i: int, depth: int) -> tuple[Match, object] | None:
    """First matching handler for outside neighbor nb on ring carbon ring_i."""
    for h in _HANDLERS:
        if h.complex and depth >= 3:
            continue
        m = h.match(mol, nb, ring_i, depth)
        if m is not None:
            return m, h
    return None


def _topology(match: Match) -> LeafTopology:
    value = match.get("z", match.get("n", 0))
    extra = match.get("o_idx", match.get("ch2", -1))
    return LeafTopology(match[ArylLeafKind], match["site"], frozenset(match["atoms"]), value,
                        frozenset(match.get("child_ring", ())), match.get("child_attach", -1),
                        match.get("child_parent", -1), extra)


def match_leaf_topology(mol: Mol, nb, ring_i: int, depth: int) -> LeafTopology | None:
    got = match_leaf(mol, nb, ring_i, depth)
    return _topology(got[0]) if got is not None else None


def match_leaf_kind(mol: Mol, nb, ring_i: int, depth: int = 1) -> ArylLeafKind | None:
    got = match_leaf(mol, nb, ring_i, depth)
    return None if got is None else got[0][ArylLeafKind]



def nested_leaf_atoms(
    mol: Mol, ring: set[int], attach: int, parent: int, depth: int = 1,
) -> set[int]:
    """Pure-topology nested atoms: ring ∪ every matched leaf, recursing into
    complex (nested ring) leaves. Topology-only; naming stays in layer3."""
    out = set(ring)
    for i in ring:
        for nb in nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            got = match_leaf(mol, nb, i, depth)
            if got is None:
                continue
            m, h = got
            out |= set(m["atoms"])
            if h.complex and m.get("child_ring"):
                out |= nested_leaf_atoms(
                    mol, m["child_ring"], m["child_attach"], m["child_parent"],
                    depth + 1,
                )
    return out



def all_outside_matched(
    mol: Mol, ring: set[int], attach: int, parent: int, depth: int,
) -> bool:
    """True if every outside (except parent link) matches a leaf; ≤3 leaves."""
    n = 0
    for i in ring:
        for nb in nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if match_leaf(mol, nb, i, depth) is None:
                return False
            n += 1
    return n <= 3
