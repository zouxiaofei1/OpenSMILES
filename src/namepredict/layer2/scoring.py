"""Composable parent-candidate scoring (simplified IUPAC P-44 seniority).

Each candidate parent dict is scored into a tuple (bigger wins):
(has_principal_fg, fg_class_rank, sides_ok, is_hetero_ring, is_carbo_ring,
 n_rings, ring_size, retained_bonus, n_unsat, n_carbons, -n_unhandled_side)

`sides_ok` sits above the ring bits on purpose: a ring candidate whose
side chains cannot be expressed as substituents (parent["n_unhandled"]>0)
must lose to a chain fallback instead of emitting a bare ring name.
"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr

# Thin frozenset wrappers derived from registry (compat / external imports).
_FG_RANK = {
    k: _kr.fg_rank(k) for k in _kr.all_kinds() if _kr.fg_rank(k)
}
_HETERO_RING = frozenset(k for k in _kr.all_kinds() if _kr.is_hetero_ring(k))
_CARBO_RING = frozenset(k for k in _kr.all_kinds() if _kr.is_carbo_ring(k))
_RETAINED = frozenset(k for k in _kr.all_kinds() if _kr.retained_bonus(k))


def _has_principal_fg(kind: str) -> int:
    return _kr.has_principal_fg(kind)


def _fg_class_rank(kind: str) -> int:
    return _kr.fg_rank(kind)


def _n_unhandled(parent: dict) -> int:
    return int(parent.get("n_unhandled") or 0)


def _sides_ok(parent: dict) -> int:
    return 0 if _n_unhandled(parent) else 1


def _is_hetero_ring(kind: str) -> int:
    return _kr.is_hetero_ring(kind)


def _is_carbo_ring(kind: str) -> int:
    return _kr.is_carbo_ring(kind)


def _n_rings(kind: str) -> int:
    return _kr.n_rings_of(kind)


def _ring_size(parent: dict, kind: str) -> int:
    if not _n_rings(kind):
        return 0
    return len(parent.get("chain") or [])


def _n_unsat(parent: dict) -> int:
    kind = parent.get("kind")
    if kind in ("polyene", "cyclopolyene"):
        return len(parent.get("double_bonds") or [])
    return 1 if kind in ("alkene", "alkyne", "cycloalkene") else 0


def _n_carbons(parent: dict) -> int:
    return int(parent.get("n_carbons") or len(parent.get("chain") or []))


def _retained_bonus(kind: str) -> int:
    return _kr.retained_bonus(kind)


def _score_parent(info: dict, parent: dict) -> tuple:
    kind = parent.get("kind") or ""
    return (
        _has_principal_fg(kind), _fg_class_rank(kind), _sides_ok(parent),
        _is_hetero_ring(kind), _is_carbo_ring(kind),
        _n_rings(kind), _ring_size(parent, kind),
        _retained_bonus(kind), _n_unsat(parent), _n_carbons(parent),
        -_n_unhandled(parent),
    )


def _better_parent(info: dict, a: dict, b: dict) -> bool:
    return _score_parent(info, a) > _score_parent(info, b)


def _pick_best(info: dict, candidates: list[dict]) -> dict | None:
    best: dict | None = None
    for cand in candidates:
        if cand is None:
            continue
        if best is None or _better_parent(info, cand, best):
            best = cand
    return best
