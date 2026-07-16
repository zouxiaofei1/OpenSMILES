"""Anthracene / 9,10-anthraquinone chain orientation (P-14 / P-25)."""
from __future__ import annotations

from namepredict.layer4.locants.adapt import ANTHRA_LOCANTS


def _anthra_loc_on(chain: list[int], attach: int) -> int:
    if attach not in chain or len(chain) != len(ANTHRA_LOCANTS):
        return 99
    loc = ANTHRA_LOCANTS[chain.index(attach)]
    return 99 if loc is None else loc


def _anthra_loc_key(chain: list[int], substituents: list) -> tuple:
    locs = sorted(_anthra_loc_on(chain, s["attach_idx"]) for s in substituents)
    return tuple(locs) if locs else ()


def _pick_anthra_chain(cands: list, key_fn) -> list[int]:
    best = cands[0]
    for cand in cands[1:]:
        if key_fn(cand) < key_fn(best):
            best = cand
    return best


def orient_anthraquinone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    cands = parent.get("anthra_chains") or [chain]
    if not cands:
        return chain
    if not substituents:
        return cands[0]
    return _pick_anthra_chain(cands, lambda c: _anthra_loc_key(c, substituents))
