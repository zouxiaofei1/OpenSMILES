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


def _prefer_lowest_locs(a, b, locs_fn, prefer_fn) -> list:
    """Pick direction with the lower locant set; ties go to prefer_fn."""
    la, lb = locs_fn(a), locs_fn(b)
    if la is None:
        return b
    if lb is None or la < lb:
        return a
    if lb < la:
        return b
    return prefer_fn(a, b)


def _prefer_lowest_bond_locs(a, b, bonds, subs, prefer_fn) -> list:
    """Pick direction with the lower bond-locant set; ties go to prefer_fn."""
    return _prefer_lowest_locs(
        a, b, lambda x: _bond_min_locs(x, bonds), lambda x, y: prefer_fn(x, y, subs),
    )


def _pair_locants(chain: list[int], cs) -> tuple[int, ...] | None:
    if not cs:
        return None
    locs = sorted(chain.index(c) + 1 for c in cs if c in chain)
    return tuple(locs) if len(locs) == len(cs) else None


def _argmin(cands, key_fn):
    """Return the candidate minimizing key_fn (first wins ties)."""
    best = cands[0]
    best_k = key_fn(best)
    for cand in cands[1:]:
        k = key_fn(cand)
        if k < best_k:
            best, best_k = cand, k
    return best


def _orient_key(chain: list[int], substituents: list) -> tuple:
    """Locant set first, then stem-alpha pairs when sets tie (P-14.5)."""
    return (tuple(_sub_locants(chain, substituents)), _stem_loc_pairs(chain, substituents))


def _prefer(a: list[int], b: list[int], subs: list) -> list[int]:
    return a if _orient_key(a, subs) <= _orient_key(b, subs) else b


def _pick_ring_by_pair_locants(chain, cs, subs, cands_fn, prefer_fn):
    """Best ring candidate by pair-locants; ties go to prefer_fn (order kept)."""
    best, best_locs = chain, _pair_locants(chain, cs)
    for cand in cands_fn(chain):
        locs = _pair_locants(cand, cs)
        if locs is None:
            continue
        if best_locs is None or locs < best_locs:
            best, best_locs = cand, locs
        elif locs == best_locs:
            best = prefer_fn(best, cand, subs)
    return best


def _table_loc_on(chain: list[int], attach: int, table) -> int:
    """Locant of `attach` by index into a fixed locant table; 99 if off-table."""
    if attach not in chain or len(chain) != len(table):
        return 99
    loc = table[chain.index(attach)]
    return 99 if loc is None else loc


def _table_loc_key(chain: list[int], substituents: list, table) -> tuple:
    locs = sorted(_table_loc_on(chain, s["attach_idx"], table) for s in substituents)
    return tuple(locs) if locs else ()


def _orient_table_cands(chain, parent, substituents, key, loc_key_fn) -> list[int]:
    """Pick among table-locant candidate chains; fall back to chain/cands[0]."""
    cands = parent.get(key) or [chain]
    if not cands:
        return chain
    if not substituents:
        return cands[0]
    return _argmin(cands, lambda c: loc_key_fn(c, substituents))


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
