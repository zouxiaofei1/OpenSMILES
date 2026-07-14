"""Leaf registry: complex handlers first, then simple."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.leaves.complex_h import COMPLEX_HANDLERS
from namepredict.layer2.leaves.protocol import Match
from namepredict.layer2.leaves.simple import SIMPLE_HANDLERS
from namepredict.layer2.leaves.topo import nb_out

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


def match_leaf_kind(mol: Mol, nb, ring_i: int, depth: int = 1) -> str | None:
    got = match_leaf(mol, nb, ring_i, depth)
    return None if got is None else got[0]["kind"]


def name_leaf(mol: Mol, m: Match, h, depth: int) -> tuple[str, str, set[int]]:
    return h.name(mol, m, depth)


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
