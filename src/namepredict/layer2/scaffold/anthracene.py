"""Retained anthracene parent (IUPAC P-25).

Linear fused three C6 aromatic rings (14 C). Unsubstituted or mono-Me/halo.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _is_linear(system: dict) -> bool:
    edges = system.get("fusion_edges") or []
    if len(edges) != 2:
        return False
    a, b = {edges[0][0], edges[0][1]}, {edges[1][0], edges[1][1]}
    return len(a & b) == 1
