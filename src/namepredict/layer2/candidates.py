"""Layer2 parent-candidate collection (scored by layer2.scoring, P-44).

FG, ring, and unsat classes each try independently (no short-circuit `or`).
Scoring then picks the best among all viable parents so `_FG_RANK`
actually arbitrates acid vs alcohol vs amine, etc.

FG parent producers are registered in `layer2.fg_producers` →
`kind_registry.fg_try_fns()`; ring parents in `layer2.ring_producers` →
`kind_registry.ring_try_fns()`; open-chain unsat hydrocarbon parents in
`layer2.unsat_producers` → `kind_registry.unsat_try_fns()` (not hand-written
tuples here).
"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr
from namepredict.layer2.parent_selector import (
    _benzene_parent,
    _longest_chain,
    _parent_dict,
)
from namepredict.layer2.ring_parent import (
    _is_benzene_core,
    _is_simple_benzene,
)
from namepredict.layer2.scoring import _pick_best


def _benzene_candidate(info: dict) -> dict | None:
    if not _is_benzene_core(info):
        return None
    cand = _benzene_parent(info)
    if not _is_simple_benzene(info):
        cand["n_unhandled"] = 1
    return cand


def _alkane_fallback(info: dict) -> dict:
    return _parent_dict(_longest_chain(info["mol"]), "alkane")


def _fg_candidates(info: dict) -> list[dict]:
    return [c for fn in _kr.fg_try_fns() if (c := fn(info)) is not None]


def _ring_candidates(info: dict) -> list[dict]:
    return [c for fn in _kr.ring_try_fns() if (c := fn(info)) is not None]


def _unsat_candidates(info: dict) -> list[dict]:
    return [c for fn in _kr.unsat_try_fns() if (c := fn(info)) is not None]


def _dedupe_parents(cands: list[dict]) -> list[dict]:
    seen: set[tuple] = set()
    out: list[dict] = []
    for c in cands:
        key = (c.get("kind"), tuple(c.get("chain") or []))
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _collect_candidates(info: dict) -> list[dict]:
    raw = (
        _fg_candidates(info) + _ring_candidates(info)
        + [_benzene_candidate(info)] + _unsat_candidates(info)
        + [_alkane_fallback(info)]
    )
    return _dedupe_parents([c for c in raw if c is not None])


def _fg_parent(info: dict) -> dict | None:
    """Compat: best FG among independent class tries."""
    return _pick_best(info, _fg_candidates(info))


def _ring_parent(info: dict) -> dict | None:
    """Compat: best ring among independent class tries."""
    return _pick_best(info, _ring_candidates(info))
