"""Layer2 parent-candidate collection (scored by layer2.scoring, P-44).

Ring parents come from `rule_driven_parent_candidates` (scaffold/core-table
driven). The former `ring_producers` / `unsat_producers` bootstrap chains in
kind_registry were removed as dead code.
"""
from __future__ import annotations

from namepredict.layer2.parent_core import _longest_chain, _parent_dict
from namepredict.layer2.parent_candidate import with_principal_group_contract
from namepredict.layer2.principal_parent import rule_driven_parent_candidates


def _alkane_fallback(info: dict) -> dict:
    return _parent_dict(_longest_chain(info["mol"]), "alkane")


def _candidate_key(candidate: dict) -> tuple:
    return candidate.get("kind"), tuple(candidate.get("chain") or [])


def _dedupe_parents(cands: list[dict]) -> list[dict]:
    seen: set[tuple] = set()
    out: list[dict] = []
    for raw in cands:
        candidate = with_principal_group_contract(raw)
        key = _candidate_key(candidate)
        if key not in seen:
            seen.add(key)
            out.append(candidate)
    return out


def _candidate_result(raw: list[dict]) -> list[dict]:
    return raw


def _collect_candidates(info: dict) -> list[dict]:
    raw = rule_driven_parent_candidates(info)
    return _dedupe_parents(_candidate_result([c for c in raw if c is not None]))
