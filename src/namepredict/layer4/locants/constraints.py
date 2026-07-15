"""P-14 soft constraint keys for ring numbering (L4 pure)."""
from __future__ import annotations

# IUPAC a-order priority rank (lower = higher priority / better low locant)
_Z_RANK: dict[int, int] = {
    8: 0,   # O
    16: 1,  # S
    34: 2,  # Se
    52: 3,  # Te
    7: 4,   # N
    15: 5,  # P
}


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


def _hetero_loc_set(
    order: tuple[int, ...] | list[int],
    hetero_atoms: list[int] | None,
) -> tuple[int, ...]:
    """Sorted locants of hetero atoms (lowest set wins)."""
    if not hetero_atoms:
        return ()
    locs = [_pos(order, a) for a in hetero_atoms]
    if any(x is None for x in locs):
        return (999,)
    return tuple(sorted(int(x) for x in locs))


def _z_of(
    atom: int,
    hetero_z: dict[int, int] | list[int] | None,
    hetero_atoms: list[int] | None,
) -> int:
    """Atomic number for atom from dict or parallel list; 999 if unknown."""
    if isinstance(hetero_z, dict):
        return int(hetero_z.get(atom, 999))
    if isinstance(hetero_z, list) and hetero_atoms:
        return _z_from_list(atom, hetero_z, hetero_atoms)
    return 999


def _z_from_list(atom: int, zs: list[int], atoms: list[int]) -> int:
    try:
        return int(zs[atoms.index(atom)])
    except (ValueError, IndexError):
        return 999


def _elem_rank_key(
    order: tuple[int, ...] | list[int],
    hetero_atoms: list[int] | None,
    hetero_z: dict[int, int] | list[int] | None,
) -> tuple[int, ...]:
    """Ranks at successive locants (O first); lower rank better at low locant."""
    if not hetero_atoms:
        return ()
    pairs = _hetero_loc_rank_pairs(order, hetero_atoms, hetero_z)
    return tuple(r for _, r in sorted(pairs))


def _hetero_loc_rank_pairs(
    order: tuple[int, ...] | list[int],
    hetero_atoms: list[int],
    hetero_z: dict[int, int] | list[int] | None,
) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for a in hetero_atoms:
        loc = _pos(order, a)
        if loc is None:
            return [(999, 999)]
        z = _z_of(a, hetero_z, hetero_atoms)
        out.append((loc, _Z_RANK.get(z, 99)))
    return out


def _multi_hetero_key(
    order: tuple[int, ...] | list[int],
    *,
    hetero_atoms: list[int] | None,
    hetero_z: dict[int, int] | list[int] | None,
    sub_attach: list[int] | None,
) -> tuple:
    """hetero set → element priority → sub set → stable."""
    return (
        _hetero_loc_set(order, hetero_atoms),
        _elem_rank_key(order, hetero_atoms, hetero_z),
        _sub_loc_set(order, sub_attach),
        _stable_tie(order),
    )


def constraint_key(
    order: tuple[int, ...] | list[int],
    mode: str,
    *,
    double_bonds: list[tuple[int, int]] | None = None,
    sub_attach: list[int] | None = None,
    hetero_atoms: list[int] | None = None,
    hetero_z: dict[int, int] | list[int] | None = None,
) -> tuple:
    """Lexicographic soft key: lower is better (P-14.4)."""
    if mode == "multi_hetero":
        return _multi_hetero_key(
            order, hetero_atoms=hetero_atoms, hetero_z=hetero_z,
            sub_attach=sub_attach,
        )
    return _carbo_or_poly_key(order, mode, double_bonds, sub_attach)


def _carbo_or_poly_key(order, mode, double_bonds, sub_attach) -> tuple:
    subs = _sub_loc_set(order, sub_attach)
    stable = _stable_tie(order)
    if mode == "poly_unsat":
        return (_ene_loc_set(order, double_bonds), subs, stable)
    return (subs, stable)


def constraints_applied(mode: str) -> tuple[str, ...]:
    if mode == "poly_unsat":
        return ("unsaturation", "substituent")
    if mode == "carbocycle_free":
        return ("substituent",)
    if mode == "multi_hetero":
        return ("hetero", "substituent")
    return ()
