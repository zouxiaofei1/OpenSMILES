"""Retained 1,3-benzothiazole / benzothiazolamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[d]thiazole. S=1, N=3 (1,3 on five-ring);
≤2 halo / methyl / CF3; mono primary amine at C2 → benzothiazolamine
(with ≤1 halo or CF3).

Core gate/subs/parent dict via fused56.Fused56Di13Spec engine.
"""
from __future__ import annotations

from namepredict.layer2.fused56 import Fused56Di13Spec, _try_di13_amine, _try_di13_fused56

_BTZ = Fused56Di13Spec(
    kind="benzothiazole",
    hetero_z=16,
    hetero_key="s_idx",
    amine_kind="benzothiazolamine",
)


def _try_benzothiazole_parent(info: dict) -> dict | None:
    return _try_di13_fused56(info, _BTZ)


def _try_benzothiazolamine_parent(info: dict) -> dict | None:
    return _try_di13_amine(info, _BTZ)
