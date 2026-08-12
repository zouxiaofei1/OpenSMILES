"""Topology-only scaffold identity shared by downstream layer facts."""
from __future__ import annotations

from dataclasses import dataclass


ScaffoldId = str


@dataclass(frozen=True)
class ScaffoldIdentity:
    id: ScaffoldId
    naming_class: str
    n_rings: int
    ring: str


def identity_of(spec) -> ScaffoldIdentity:
    return ScaffoldIdentity(spec.id, spec.naming_class, spec.n_rings, spec.ring)
