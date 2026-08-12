"""Anthracene / 9,10-anthraquinone chain orientation (P-14 / P-25)."""
from __future__ import annotations

from namepredict.layer4._chain_orient import _orient_table_cands, _table_loc_key
from namepredict.layer4.locants.adapt import ANTHRA_LOCANTS


def orient_anthraquinone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_table_cands(
        chain, parent, substituents, "anthra_chains",
        lambda c, s: _table_loc_key(c, s, ANTHRA_LOCANTS),
    )
