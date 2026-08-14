"""多烯位次，供 layer4 FG 位次计算使用（P-31.1）。"""
from __future__ import annotations

from namepredict.layer4._chain_orient import _bond_min_locs


def ene_locants(oriented: dict) -> list[int] | None:
    bonds = oriented.get("double_bonds")
    if not bonds:
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None
