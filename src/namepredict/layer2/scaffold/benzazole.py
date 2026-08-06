"""Retained 1,3-benzazole / benzazolamine parents (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[d]thiazole (S=1, N=3) and benzo[d]oxazole
(O=1, N=3).  Engine via fused56; Di13Spec + ScaffoldSpec from
scaffold.builders.fused56.
"""
from __future__ import annotations

from namepredict.layer2.scaffold.fused56 import _try_di13_amine, _try_di13_fused56
from namepredict.layer2.scaffold.builders.fused56 import BOX_DI13, BTZ_DI13


def _try_benzothiazole_parent(info: dict) -> dict | None:
    return _try_di13_fused56(info, BTZ_DI13)


def _try_benzothiazolamine_parent(info: dict) -> dict | None:
    return _try_di13_amine(info, BTZ_DI13)


def _try_benzoxazole_parent(info: dict) -> dict | None:
    return _try_di13_fused56(info, BOX_DI13)


def _try_benzoxazolamine_parent(info: dict) -> dict | None:
    return _try_di13_amine(info, BOX_DI13)
