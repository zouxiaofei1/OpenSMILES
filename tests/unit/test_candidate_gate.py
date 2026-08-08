# Layer: L2
"""Candidate-gate contract: producer scope must never leak globally."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.arene_carbonyl import benzene_polycarboxylic_gate
from namepredict.layer2.candidate_gate import (
    CandidateGate, GateScope, GateStatus, gate_result, global_reject, pass_gate,
    scoped_reject,
)
from namepredict.layer2.candidates import _collect_candidates, _polyacid_gates


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _kinds(smiles: str) -> set[str]:
    return {candidate["kind"] for candidate in _collect_candidates(_info(smiles))}


@pytest.mark.parametrize("smiles", [
    "O=C(O)CC(Cc1ccccc1)C(=O)O",
    "O=C(O)C1CCCCC1C(=O)O",
    "O=C(O)c1ccncc1C(=O)O",
])
def test_benzene_gate_passes_non_scope_structures(smiles: str) -> None:
    gate = benzene_polycarboxylic_gate(_info(smiles))
    assert gate.status is GateStatus.PASS
    assert gate.scope is GateScope.BENZENE_POLYCARBOXYLIC


def test_benzene_gate_rejects_unsupported_structure_in_its_scope() -> None:
    gate = benzene_polycarboxylic_gate(_info("O=C(O)c1ccccc1C(=O)OC"))
    assert gate.status is GateStatus.SCOPED_REJECT
    assert gate.reason == "unsupported_benzene_polyacid"
    assert _kinds("O=C(O)c1ccccc1C(=O)OC") == {"unsupported_polycarboxylic"}


@pytest.mark.parametrize("smiles", [
    "O=C(O)CC(Cc1ccccc1)C(=O)O",
    "O=C(O)C1CCCCC1C(=O)O",
    "O=C(O)c1ccncc1C(=O)O",
])
def test_benzene_gate_does_not_terminate_other_scopes(smiles: str) -> None:
    assert all(gate.status is GateStatus.PASS for gate in _polyacid_gates(_info(smiles)))
    assert "unsupported_polycarboxylic" not in _kinds(smiles)


def test_pass_gate_keeps_candidates() -> None:
    candidates = [{"kind": "acid"}]
    assert gate_result([pass_gate(GateScope.BENZENE_POLYCARBOXYLIC)], candidates) == (candidates, None)


def test_scoped_reject_keeps_unowned_candidate() -> None:
    scope = GateScope.BENZENE_POLYCARBOXYLIC
    candidates = [{"kind": "other", "gate_principal": True}, {"kind": "benzene", "gate_dependencies": (scope,)}]
    assert gate_result([scoped_reject(scope, "bad")], candidates) == ([candidates[0]], None)


def test_two_scoped_rejections_remove_only_their_owned_candidates() -> None:
    benzene = GateScope.BENZENE_POLYCARBOXYLIC
    chain = GateScope.OPEN_CHAIN_POLYCARBOXYLIC
    candidates = [
        {"gate_dependencies": (benzene,)}, {"gate_dependencies": (chain,)},
        {"kind": "other", "gate_principal": True},
    ]
    gates = [scoped_reject(benzene, "benzene"), scoped_reject(chain, "chain")]
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
    gate = global_reject(GateScope.BENZENE_POLYCARBOXYLIC, "global")
    assert gate_result([gate], [candidate]) == ([], "global")


def test_gate_scope_is_closed_enum() -> None:
    with pytest.raises(TypeError):
        CandidateGate(GateStatus.PASS, "unknown_scope")
