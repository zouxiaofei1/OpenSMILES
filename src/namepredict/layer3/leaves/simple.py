"""Simple (non-recursive) leaf match handlers: halo/Me/alkoxy/nitro/OH/NH2/CF3.

Match-only registry entries; naming lives in L3 (name_as_substituent).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer3.leaves import topo
from namepredict.layer3.leaves.protocol import Match


class _FnHandler:
    """Adapter: match callable + complex flag."""

    def __init__(self, match_fn, complex: bool = False):
        self.complex = complex
        self._match = match_fn

    def match(self, mol: Mol, nb, ring_i: int, depth: int) -> Match | None:
        return self._match(mol, nb, ring_i, depth)


SIMPLE_HANDLERS = [
    _FnHandler(topo.match_halo),
    _FnHandler(topo.match_amino),
    _FnHandler(topo.match_cyano),
    _FnHandler(topo.match_alkoxy),
    _FnHandler(topo.match_hydroxy),
    _FnHandler(topo.match_me),
    _FnHandler(topo.match_n_alkyl),
    _FnHandler(topo.match_methylthio),
    _FnHandler(topo.match_nitro),
    _FnHandler(topo.match_cf3),
]
