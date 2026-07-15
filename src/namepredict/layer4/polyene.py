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


def orient_alkenedioic(chain, parent, subs, prefer_fn, mono_fn):
    """Multi-ene: lowest ene set; mono C=C: reuse mono alkene orienter."""
    if parent.get("double_bonds"):
        return orient_polyene(chain, parent, subs, prefer_fn)
    return mono_fn(chain, parent, subs)


def ene_locants(oriented: dict) -> list[int] | None:
    bonds = oriented.get("double_bonds")
    if not bonds:
        return None
    kind = oriented.get("kind")
    if kind not in ("polyene", "alkenol", "alkenoic_acid", "alkenedioic"):
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None


def _pick_oh_orient(a: list[int], b: list[int], oh: int) -> list[int] | None:
    la, lb = a.index(oh) + 1, b.index(oh) + 1
    if la < lb:
        return a
    return b if lb < la else None


def orient_alkenol(chain: list[int], parent: dict, subs: list, prefer_fn) -> list[int]:
    """OH lowest first; multi-ene uses polyene set as secondary (P-31.1)."""
    oh = parent.get("oh_c_idx")
    if oh is None or oh not in chain:
        return orient_polyene(chain, parent, subs, prefer_fn)
    picked = _pick_oh_orient(list(chain), list(reversed(chain)), oh)
    return picked if picked is not None else orient_polyene(chain, parent, subs, prefer_fn)
