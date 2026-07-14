"""LeafHandler protocol types for aryl ring outside-neighbors."""
from __future__ import annotations

from typing import Callable, Protocol

from rdkit.Chem import Mol

# Match: topology result for one outside neighbor on a ring carbon.
# keys: kind, atoms, site (ring carbon), optional child fields for complex leaves
Match = dict

# ParentCtx for match/name
# mol, ring, ring_i (attachment carbon on parent ring), depth, max_depth
Ctx = dict


class LeafHandler(Protocol):
    key: str
    alpha_key: str
    complex: bool

    def match(self, mol: Mol, nb, ring_i: int, depth: int) -> Match | None: ...

    def name(self, mol: Mol, m: Match, depth: int) -> tuple[str, str, set[int]]: ...


def make_match(kind: str, atoms: set[int], site: int, **extra) -> Match:
    return {"kind": kind, "atoms": atoms, "site": site, **extra}
