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


def _fg_pos(chain: list[int], c: int | None) -> int | None:
    return None if c is None or c not in chain else chain.index(c) + 1


def prefer_unsat_if_fg_tie(base, parent, subs, fg_key, prefer_ene_fn):
    """When FG locant ties both ways, prefer lower ene/yne locant (P-31.1)."""
    ends = parent.get("double_bond") or parent.get("triple_bond")
    if not ends:
        return base
    rev = list(reversed(base))
    if _fg_pos(base, parent.get(fg_key)) != _fg_pos(rev, parent.get(fg_key)):
        return base
    return prefer_ene_fn(base, rev, ends, subs)


def orient_cyclopolyene(chain: list[int], parent: dict, subs: list) -> list[int]:
    """Ring polyene: lowest endocyclic ene set via locant engine (poly_unsat)."""
    from namepredict.layer4.locants.engine import choose_numbering
    bonds = parent.get("double_bonds") or []
    attach = [s.get("attach_idx") for s in subs if s.get("attach_idx") is not None]
    plan = choose_numbering(
        chain, "poly_unsat", double_bonds=bonds, sub_attach=attach,
        scaffold_id="cyclopolyene",
    )
    return list(plan.atom_order) if plan.atom_order else chain


def _rotate_to(chain: list[int], atom: int) -> list[int]:
    if atom not in chain:
        return chain
    i = chain.index(atom)
    return chain[i:] + chain[:i]


def _ene_base_from(chain: list[int], a: int, b: int) -> list[int] | None:
    base = _rotate_to(chain, a)
    if len(base) > 1 and base[1] != b:
        base = _rotate_to(list(reversed(chain)), a)
    return base if len(base) > 1 and base[1] == b else None


def _cyclo_ene_bases(chain: list[int], ends: tuple) -> list[list[int]]:
    """Ring bases with C=C as 1–2, starting from either end."""
    out: list[list[int]] = []
    for a, b in (ends, (ends[1], ends[0])):
        base = _ene_base_from(chain, a, b) if a in chain else None
        if base is not None:
            out.append(base)
    return out


def orient_cycloalkene(chain: list[int], parent: dict, subs: list, prefer_fn) -> list[int]:
    """Double bond at 1–2; choose direction by lowest substituent set."""
    ends = parent.get("double_bond")
    if not ends or ends[0] not in chain:
        return chain
    bases = _cyclo_ene_bases(chain, ends)
    best = bases[0] if bases else chain
    for cand in bases[1:]:
        best = prefer_fn(best, cand, subs)
    return best


def orient_ring_fg_ene(
    chain: list[int], parent: dict, subs: list, key: str, prefer_fn, prefer_ene_fn,
) -> list[int]:
    """FG@1 both directions; prefer lower ene locant, then sub set."""
    c = parent.get(key)
    if c is None or c not in chain:
        return chain
    base, rev = _rotate_to(chain, c), _rotate_to(list(reversed(chain)), c)
    ends = parent.get("double_bond")
    if ends:
        return prefer_ene_fn(base, rev, ends, subs)
    return prefer_fn(base, rev, subs)
