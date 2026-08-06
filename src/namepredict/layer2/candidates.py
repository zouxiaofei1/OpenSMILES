"""Layer2 parent-candidate collection (scored by layer2.scoring, P-44).

Ring parents in `layer2.ring_producers` → `kind_registry.ring_try_fns()`;
open-chain unsat hydrocarbon parents in `layer2.unsat_producers` →
`kind_registry.unsat_try_fns()`.
"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr
from namepredict.layer2.parent_core import _longest_chain, _parent_dict
from namepredict.layer2.parent_candidate import with_principal_group_contract
from namepredict.layer2.principal_parent import rule_driven_parent_candidates
from namepredict.layer2.arene_carbonyl import benzene_polycarboxylic_gate
from namepredict.layer2.candidate_gate import CandidateGate, GateScope, gate_result
from namepredict.layer2.carboxymethyl_diacid import is_carboxymethyl_diacid
from namepredict.layer2.polycarboxylic import polycarboxylic_gate
from namepredict.layer2.cyclo_polycarboxylic import cycloalkane_polycarboxylic_gate


def _alkane_fallback(info: dict) -> dict:
    return _parent_dict(_longest_chain(info["mol"]), "alkane")


def _unsat_candidates(info: dict) -> list[dict]:
    return [c for fn in _kr.unsat_try_fns() if (c := fn(info)) is not None]


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


def _unsupported_polyacid(reason: str | None) -> list[dict]:
    return [{"kind": "unsupported_polycarboxylic", "chain": [], "n_carbons": 0,
             "unsupported_reason": reason}]


def _polyacid_gates(info: dict) -> list[CandidateGate]:
    gates = [benzene_polycarboxylic_gate(info), polycarboxylic_gate(info), cycloalkane_polycarboxylic_gate(info)]
    return gates if not is_carboxymethyl_diacid(info) else gates[:1]


_CANDIDATE_POLICIES = {
    "acid": ((GateScope.OPEN_CHAIN_POLYCARBOXYLIC, GateScope.BENZENE_POLYCARBOXYLIC, GateScope.CYCLOALKANE_POLYCARBOXYLIC), False),
    "alkane": ((GateScope.OPEN_CHAIN_POLYCARBOXYLIC, GateScope.BENZENE_POLYCARBOXYLIC, GateScope.CYCLOALKANE_POLYCARBOXYLIC), False),
    "benzene": ((GateScope.BENZENE_POLYCARBOXYLIC,), False),
    "benzene_polycarboxylic": ((GateScope.BENZENE_POLYCARBOXYLIC,), True),
    "cycloalkane_polycarboxylic": ((GateScope.CYCLOALKANE_POLYCARBOXYLIC,), True),
    "benzoic": ((), True),
    "diacid": ((GateScope.BENZENE_POLYCARBOXYLIC,), False),
    "polycarboxylic": ((GateScope.OPEN_CHAIN_POLYCARBOXYLIC,), True),
}


def _policy_fields(candidate: dict) -> dict:
    dependencies, principal = _CANDIDATE_POLICIES.get(candidate.get("kind"), ((), False))
    return {"gate_dependencies": dependencies, "gate_principal": principal}


def _owned_candidate(candidate: dict) -> dict:
    return {**candidate, **_policy_fields(candidate)}


def _candidate_result(info: dict, raw: list[dict]) -> list[dict]:
    kept, reason = gate_result(_polyacid_gates(info), [_owned_candidate(c) for c in raw])
    return _unsupported_polyacid(reason) if reason is not None else kept


def _principal_candidates(info: dict) -> list[dict]:
    return rule_driven_parent_candidates(info)


def _collect_candidates(info: dict) -> list[dict]:
    raw = _principal_candidates(info)
    return _dedupe_parents(_candidate_result(info, [c for c in raw if c is not None]))
