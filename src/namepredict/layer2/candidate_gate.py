"""Typed Layer 2 candidate-gating contract."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class GateStatus(Enum):
    PASS = "pass"
    SCOPED_REJECT = "scoped_reject"
    GLOBAL_REJECT = "global_reject"


class GateScope(Enum):
    OPEN_CHAIN_POLYCARBOXYLIC = "open_chain_polycarboxylic"


@dataclass(frozen=True)
class CandidateGate:
    """A producer decision constrained to a registered candidate scope."""

    status: GateStatus
    scope: GateScope
    reason: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.scope, GateScope):
            raise TypeError("scope must be a GateScope")


def pass_gate(scope: GateScope) -> CandidateGate:
    return CandidateGate(GateStatus.PASS, scope)


def scoped_reject(scope: GateScope, reason: str) -> CandidateGate:
    return CandidateGate(GateStatus.SCOPED_REJECT, scope, reason)


def global_reject(scope: GateScope, reason: str) -> CandidateGate:
    return CandidateGate(GateStatus.GLOBAL_REJECT, scope, reason)


def gate_result(gates: list[CandidateGate], candidates: list[dict]) -> tuple[list[dict], str | None]:
    """Interpret claimed scopes, candidate dependencies, and principals separately."""
    global_gate = next((gate for gate in gates if gate.status is GateStatus.GLOBAL_REJECT), None)
    if global_gate is not None:
        return [], global_gate.reason
    return _scoped_result(gates, candidates)


def _scoped_result(gates: list[CandidateGate], candidates: list[dict]) -> tuple[list[dict], str | None]:
    rejected = [gate for gate in gates if gate.status is GateStatus.SCOPED_REJECT]
    kept = [candidate for candidate in candidates if not _depends_on(candidate, rejected)]
    if not rejected or _has_independent_principal(kept):
        return kept, None
    return [], rejected[0].reason


def _depends_on(candidate: dict, gates: list[CandidateGate]) -> bool:
    dependencies = candidate.get("gate_dependencies", ())
    return any(gate.scope in dependencies for gate in gates)


def _has_independent_principal(candidates: list[dict]) -> bool:
    return any(candidate.get("gate_principal") for candidate in candidates)
