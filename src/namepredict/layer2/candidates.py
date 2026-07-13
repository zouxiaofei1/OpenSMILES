"""Layer2 parent-candidate collection (scored by layer2.scoring, P-44).

Keeps `parent_selector` waterfall try-functions intact; each path now
contributes a candidate instead of short-circuiting. A benzene core is
always a candidate — when its sides are not fully expressible as
substituents it carries n_unhandled so scoring prefers a chain fallback.
"""
from __future__ import annotations

from namepredict.layer2.parent_selector import (
    _benzene_parent,
    _fg_parent,
    _longest_chain,
    _parent_dict,
    _ring_parent,
    _unsat_parent,
)
from namepredict.layer2.ring_parent import _is_benzene_core, _is_simple_benzene


def _benzene_candidate(info: dict) -> dict | None:
    if not _is_benzene_core(info):
        return None
    cand = _benzene_parent(info)
    if not _is_simple_benzene(info):
        cand["n_unhandled"] = 1
    return cand


def _alkane_fallback(info: dict) -> dict:
    return _parent_dict(_longest_chain(info["mol"]), "alkane")


def _collect_candidates(info: dict) -> list[dict]:
    cands = [
        _fg_parent(info), _ring_parent(info), _benzene_candidate(info),
        _unsat_parent(info), _alkane_fallback(info),
    ]
    return [c for c in cands if c is not None]
