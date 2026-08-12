"""Shared chain/ring orientation primitives for layer4 orienters.

Extracted verbatim from layer4/numbering.py, layer4/polyene.py and
layer4/sat_hetero_orient.py, which duplicated these helpers.  Only pure
helpers with identical behaviour were merged; ring-candidate generators whose
iteration order affects tie-breaking are kept at their call sites.
"""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key


def _chain_pos(chain: list[int], c: int | None) -> int | None:
    return None if c is None or c not in chain else chain.index(c) + 1


def _edge_locants(chain: list[int], edge) -> tuple[int, int] | None:
    if not edge or edge[0] not in chain or edge[1] not in chain:
        return None
    return chain.index(edge[0]) + 1, chain.index(edge[1]) + 1


def _edge_min_locant(chain: list[int], edge) -> int | None:
    ends = _edge_locants(chain, edge)
    return min(ends) if ends else None


def _bond_min_locs(chain: list[int], bonds) -> tuple[int, ...] | None:
    if not bonds:
        return None
    locs = [_edge_min_locant(chain, b) for b in bonds]
    if any(x is None for x in locs):
        return None
    return tuple(sorted(int(x) for x in locs))


def _prefer_lowest_bond_locs(a, b, bonds, subs, prefer_fn) -> list:
    """Pick direction with the lower bond-locant set; ties go to prefer_fn."""
    la, lb = _bond_min_locs(a, bonds), _bond_min_locs(b, bonds)
    if la is None:
        return b
    if lb is None or la < lb:
        return a
    if lb < la:
        return b
    return prefer_fn(a, b, subs)


def _pair_locants(chain: list[int], cs) -> tuple[int, ...] | None:
    if not cs:
        return None
    locs = sorted(chain.index(c) + 1 for c in cs if c in chain)
    return tuple(locs) if len(locs) == len(cs) else None


def _rotate_to(chain: list[int], atom: int) -> list[int]:
    if atom not in chain:
        return chain
    i = chain.index(atom)
    return chain[i:] + chain[:i]


def _sub_locants(chain: list[int], substituents: list) -> list[int]:
    # ester O-side alkyl attaches to the ester O (not on chain) and carries no
    # parent locant — exclude it from orientation comparison.
    return sorted(
        chain.index(s["attach_idx"]) + 1
        for s in substituents if s["attach_idx"] in chain
    )


def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    return sorted(
        (alkyl_alpha_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents if s["attach_idx"] in chain
    )
