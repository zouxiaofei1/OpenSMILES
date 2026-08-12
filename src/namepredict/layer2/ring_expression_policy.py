"""Capability policy for typed principal expressions on ring scaffolds."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG
from namepredict.layer2.identity import ScaffoldIdentity


@dataclass(frozen=True)
class RingExpressionPolicy:
    naming_classes: frozenset[str]
    group_class: FG
    relations: frozenset[str]
    max_multiplicity: int


_POLICIES = (
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ALCOHOL, frozenset({"in_skeleton"}), 3),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.KETONE, frozenset({"in_skeleton"}), 2),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.AMINE, frozenset({"in_skeleton"}), 1),
    RingExpressionPolicy(frozenset({"carbocycle"}), FG.ACID, frozenset({"exocyclic"}), 10),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.ALCOHOL, frozenset({"in_skeleton"}), 2),
    RingExpressionPolicy(frozenset({"mono_carbo"}), FG.AMINE, frozenset({"in_skeleton"}), 1),
    RingExpressionPolicy(frozenset({"naph_family"}), FG.ALCOHOL, frozenset({"in_skeleton"}), 1),
)


def supports_ring_expression(scaffold: ScaffoldIdentity, facts) -> bool:
    return any(scaffold.naming_class in policy.naming_classes
               and facts.group_class is policy.group_class
               and facts.relation.value in policy.relations
               and facts.multiplicity <= policy.max_multiplicity for policy in _POLICIES)
