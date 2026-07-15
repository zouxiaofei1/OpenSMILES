"""P-14 soft constraint keys for ring numbering (L4 pure)."""
from __future__ import annotations


def _pos(order: tuple[int, ...] | list[int], atom: int) -> int | None:
    try:
        return list(order).index(atom) + 1
    except ValueError:
        return None


def _sub_loc_set(
    order: tuple[int, ...] | list[int],
    sub_attach: list[int] | None,
) -> tuple[int, ...]:
    if not sub_attach:
        return ()
    locs = [_pos(order, a) for a in sub_attach]
    if any(x is None for x in locs):
        return (999,)
    return tuple(sorted(int(x) for x in locs))


def _adj_loc(lo: int, hi: int, n: int) -> int:
    """Consecutive → lo; wrap-around 1–n → n; else invalid."""
    if hi - lo == 1:
        return lo
    return n if lo == 1 and hi == n else 999


def _bond_loc(order: list[int], a: int, b: int) -> int | None:
    """Lower locant if consecutive; wrap-around uses n (disfavored)."""
    pa, pb = _pos(order, a), _pos(order, b)
    if pa is None or pb is None:
        return None
    return _adj_loc(min(pa, pb), max(pa, pb), len(order))


def _ene_loc_set(
    order: tuple[int, ...] | list[int],
    double_bonds: list[tuple[int, int]] | None,
) -> tuple[int, ...]:
    """Lower endpoint locant per endocyclic C=C, sorted (P-31.1)."""
    if not double_bonds:
        return ()
    seq = list(order)
    locs = [_bond_loc(seq, a, b) for a, b in double_bonds]
    if any(x is None for x in locs):
        return (999,)
    return tuple(sorted(int(x) for x in locs))


def _stable_tie(order: tuple[int, ...] | list[int]) -> tuple:
    """Prefer original direction (smaller first atom, then tuple order)."""
    o = tuple(order)
    return (o[0] if o else 0, o)


def constraint_key(
    order: tuple[int, ...] | list[int],
    mode: str,
    *,
    double_bonds: list[tuple[int, int]] | None = None,
    sub_attach: list[int] | None = None,
) -> tuple:
    """Lexicographic soft key: lower is better (P-14.4)."""
    subs = _sub_loc_set(order, sub_attach)
    stable = _stable_tie(order)
    if mode == "poly_unsat":
        return (_ene_loc_set(order, double_bonds), subs, stable)
    # carbocycle_free and default: substituent set primary
    return (subs, stable)


def constraints_applied(mode: str) -> tuple[str, ...]:
    if mode == "poly_unsat":
        return ("unsaturation", "substituent")
    if mode == "carbocycle_free":
        return ("substituent",)
    return ()
