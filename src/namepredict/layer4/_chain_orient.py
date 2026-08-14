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


def _pair_locants(chain: list[int], cs) -> tuple[int, ...] | None:
    if not cs:
        return None
    locs = sorted(chain.index(c) + 1 for c in cs if c in chain)
    return tuple(locs) if len(locs) == len(cs) else None


def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    return sorted(
        (alkyl_alpha_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents if s["attach_idx"] in chain
    )
