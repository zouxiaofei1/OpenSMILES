"""Simple (non-recursive) leaf handlers: halo/Me/alkoxy/nitro/OH/NH2/CF3."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import _ALKOXY_EN, _ALKOXY_ZH
from namepredict.layer2.leaves import topo
from namepredict.layer2.leaves.protocol import Match

_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}


class _FnHandler:
    """Adapter: match/name callables + metadata."""

    def __init__(self, key: str, alpha: str, match_fn, name_fn, complex: bool = False):
        self.key = key
        self.alpha_key = alpha
        self.complex = complex
        self._match = match_fn
        self._name = name_fn

    def match(self, mol: Mol, nb, ring_i: int, depth: int) -> Match | None:
        return self._match(mol, nb, ring_i, depth)

    def name(self, mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
        return self._name(mol, m, depth)


def _name_halo(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    z = m["z"]
    return _HALO_EN[z], _HALO_ZH[z], set(m["atoms"])


def _name_me(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    return "methyl", "甲基", set(m["atoms"])


def _name_cf3(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    return "(trifluoromethyl)", "三氟甲基", set(m["atoms"])


def _name_nitro(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    return "nitro", "硝基", set(m["atoms"])


def _name_hydroxy(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    return "hydroxy", "羟基", set(m["atoms"])


def _name_amino(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    return "amino", "氨基", set(m["atoms"])


def _name_alkoxy(mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]:
    n = m["n"]
    return _ALKOXY_EN[n], _ALKOXY_ZH[n], set(m["atoms"])


SIMPLE_HANDLERS = [
    _FnHandler("halo", "halo", topo.match_halo, _name_halo),
    _FnHandler("amino", "amino", topo.match_amino, _name_amino),
    _FnHandler("alkoxy", "alkoxy", topo.match_alkoxy, _name_alkoxy),
    _FnHandler("hydroxy", "hydroxy", topo.match_hydroxy, _name_hydroxy),
    _FnHandler("me", "methyl", topo.match_me, _name_me),
    _FnHandler("nitro", "nitro", topo.match_nitro, _name_nitro),
    _FnHandler("cf3", "trifluoromethyl", topo.match_cf3, _name_cf3),
]
