"""Retained 1,3-benzoxazole / benzoxazolamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[d]oxazole. O=1, N=3 (1,3 on five-ring);
≤2 halo / methyl / CF3; mono primary amine at C2 → benzoxazolamine
(with ≤1 halo or CF3).

Core gate/subs/parent dict via fused56.Fused56Di13Spec engine.
"""
from __future__ import annotations

from namepredict.layer2.fused56 import Fused56Di13Spec, _try_di13_amine, _try_di13_fused56

_BOX = Fused56Di13Spec(
    kind="benzoxazole",
    hetero_z=8,
    hetero_key="o_idx",
    amine_kind="benzoxazolamine",
)


def _try_benzoxazole_parent(info: dict) -> dict | None:
    return _try_di13_fused56(info, _BOX)


def _try_benzoxazolamine_parent(info: dict) -> dict | None:
    return _try_di13_amine(info, _BOX)
