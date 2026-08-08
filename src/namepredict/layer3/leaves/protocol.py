"""Aryl leaf-kind enum (leaf matching retired in 446ce69; enum kept for typing)."""
from __future__ import annotations

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
