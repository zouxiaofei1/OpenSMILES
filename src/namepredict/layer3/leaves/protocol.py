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


def make_match(kind: ArylLeafKind, atoms: set[int], site: int, **extra) -> Match:
    return {ArylLeafKind: kind, "atoms": atoms, "site": site, **extra}
