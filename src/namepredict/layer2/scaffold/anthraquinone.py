"""9,10-Anthraquinone retained parent (IUPAC P-25 / P-64).

Linear fused three C6 rings (14 C) with exactly two meso ring ketones at 9,10.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.anthracene import _is_linear
from namepredict.layer2.scaffold.ring_parent import (
    _dbl_o_idx,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _aq_shape_ok(s: dict) -> bool:
    return (
        s.get("n_rings") == 3
        and s.get("n_atoms") == 14
        and s.get("topology") == "fused"
        and not s.get("hetero_atoms")
        and len(s.get("fusion_edges") or []) == 2
    )


def _aq_system(info: dict) -> dict | None:
    for s in info.get("ring_systems") or []:
        if _aq_shape_ok(s) and _is_linear(s):
            return s
    return None
