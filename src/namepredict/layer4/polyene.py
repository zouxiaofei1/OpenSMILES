"""Polyene orientation and multi-ene locants (P-31.1)."""
from __future__ import annotations


def _pair_ends(chain: list[int], pair) -> tuple[int, int] | None:
    if not pair or pair[0] not in chain or pair[1] not in chain:
        return None
    return chain.index(pair[0]) + 1, chain.index(pair[1]) + 1


def _min_loc(chain: list[int], pair) -> int | None:
    ends = _pair_ends(chain, pair)
    return min(ends) if ends else None


def _bond_min_locs(chain: list[int], bonds) -> tuple[int, ...] | None:
    if not bonds:
        return None
    locs = [_min_loc(chain, b) for b in bonds]
    if any(x is None for x in locs):
        return None
    return tuple(sorted(int(x) for x in locs))


def _prefer(a: list[int], b: list[int], bonds, subs: list, prefer_fn) -> list[int]:
    la, lb = _bond_min_locs(a, bonds), _bond_min_locs(b, bonds)
    if la is None:
        return b
    if lb is None or la < lb:
        return a
    if lb < la:
        return b
    return prefer_fn(a, b, subs)


def orient_polyene(chain: list[int], parent: dict, subs: list, prefer_fn) -> list[int]:
    bonds = parent.get("double_bonds") or []
    if not bonds:
        return chain
    return _prefer(chain, list(reversed(chain)), bonds, subs, prefer_fn)


def ene_locants(oriented: dict) -> list[int] | None:
    if oriented.get("kind") != "polyene":
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], oriented.get("double_bonds"))
    return list(locs) if locs else None
