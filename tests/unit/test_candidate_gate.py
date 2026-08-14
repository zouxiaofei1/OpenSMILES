# Layer: L2
"""Candidate-gate contract: producer scope must never leak globally."""
from __future__ import annotations

import pytest

from namepredict.layer2.candidate_gate import (
    CandidateGate, GateScope, GateStatus, gate_result, global_reject, pass_gate,
    scoped_reject,
)


def test_pass_gate_keeps_candidates() -> None:
    candidates = [{"kind": "acid"}]
    assert gate_result([pass_gate(GateScope.OPEN_CHAIN_POLYCARBOXYLIC)], candidates) == (candidates, None)


def test_scoped_reject_keeps_unowned_candidate() -> None:
    scope = GateScope.OPEN_CHAIN_POLYCARBOXYLIC
    candidates = [{"kind": "other", "gate_principal": True}, {"kind": "pyridine", "gate_dependencies": (scope,)}]
    assert gate_result([scoped_reject(scope, "bad")], candidates) == ([candidates[0]], None)


def test_two_scoped_rejections_remove_only_their_owned_candidates() -> None:
    scope = GateScope.OPEN_CHAIN_POLYCARBOXYLIC
    candidates = [
        {"gate_dependencies": (scope,)}, {"gate_dependencies": (scope,)},
        {"kind": "other", "gate_principal": True},
    ]
    gates = [scoped_reject(scope, "a"), scoped_reject(scope, "b")]
    assert gate_result(gates, candidates) == ([candidates[2]], None)


def test_scoped_reject_rejects_unowned_generic_fallback() -> None:
    candidate = {"kind": "fallback"}
    gate = scoped_reject(GateScope.OPEN_CHAIN_POLYCARBOXYLIC, "bad")
    assert gate_result([gate], [candidate]) == ([], "bad")


def test_scoped_reject_terminates_only_when_no_alternative_remains() -> None:
    scope = GateScope.OPEN_CHAIN_POLYCARBOXYLIC
    candidates = [{"kind": "polyacid", "gate_scopes": (scope,)}]
    assert gate_result([scoped_reject(scope, "bad")], candidates) == ([], "bad")


def test_global_reject_terminates_independent_of_candidate_scope() -> None:
    candidate = {"kind": "other"}
    gate = global_reject(GateScope.OPEN_CHAIN_POLYCARBOXYLIC, "global")
    assert gate_result([gate], [candidate]) == ([], "global")


def test_gate_scope_is_closed_enum() -> None:
    with pytest.raises(TypeError):
        CandidateGate(GateStatus.PASS, "unknown_scope")
