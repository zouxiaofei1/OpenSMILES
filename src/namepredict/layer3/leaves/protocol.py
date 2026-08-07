"""LeafHandler protocol types for aryl ring outside-neighbors."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

class ArylLeafKind(Enum):
    HALOGEN = auto()
    AMINO = auto()
    CYANO = auto()
    ALKOXY = auto()
    HYDROXY = auto()
    METHYL = auto()
    N_ALKYL = auto()
    METHYLSULFANYL = auto()
    NITRO = auto()
    TRIFLUOROMETHYL = auto()
    PHENYL = auto()
    PHENOXY = auto()
    BENZYL = auto()


@dataclass(frozen=True)
class LeafTopology:
    kind: ArylLeafKind
    site: int
    atoms: frozenset[int]
    value: int
    child_ring: frozenset[int]
    child_attach: int
    child_parent: int
    extra_atom: int


# Match: topology result for one outside neighbor on a ring carbon.
# keys: kind, atoms, site (ring carbon), optional child fields for complex leaves
Match = dict


def make_match(
    kind: ArylLeafKind, atoms: set[int], site: int,
    *, z: int | None = None, n: int | None = None, o_idx: int | None = None,
    ch2: int | None = None, child_ring: set[int] | None = None,
    child_attach: int | None = None, child_parent: int | None = None,
) -> Match:
    m: Match = {ArylLeafKind: kind, "atoms": atoms, "site": site}
    if z is not None:
        m["z"] = z
    if n is not None:
        m["n"] = n
    if o_idx is not None:
        m["o_idx"] = o_idx
    if ch2 is not None:
        m["ch2"] = ch2
    if child_ring is not None:
        m["child_ring"] = child_ring
    if child_attach is not None:
        m["child_attach"] = child_attach
    if child_parent is not None:
        m["child_parent"] = child_parent
    return m
