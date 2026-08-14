"""仅含拓扑的骨架同一性，供下游各层 facts 共享。"""
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
