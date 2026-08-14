"""Adapt parent-owned numbering facts into an L4 NumberingPlan.

L4 accepts only ordinary parent dictionaries and never imports L2.  The
parent selector decides whether a scaffold has a plan and supplies aligned
label facts; absent facts deliberately mean no retained-scaffold plan.
"""
from __future__ import annotations

from collections.abc import Mapping

from namepredict.layer4.locants.plan import NumberingPlan, locant, make_plan


def _fact_fields(facts: object) -> tuple[str, tuple[str, ...]] | None:
    if not isinstance(facts, Mapping):
        return None
    scaffold_id, labels = facts.get("scaffold_id"), facts.get("labels")
    if not isinstance(scaffold_id, str) or not isinstance(labels, (tuple, list)):
        return None
    labels = tuple(labels)
    return (scaffold_id, labels) if all(isinstance(label, str) for label in labels) else None


def plan_from_parent(chain: list[int], parent: Mapping) -> NumberingPlan | None:
    """Build a plan from L2 materialized facts without interpreting scaffold kind."""
    fields = _fact_fields(parent.get("numbering_scaffold"))
    if fields is None or len(chain) != len(fields[1]):
        return None
    return make_plan(fields[0], tuple(chain), fields[1])


def plan_from_chain(
    chain: list[int], kind: str | None, facts: Mapping | None = None, *, required: bool = False,
) -> NumberingPlan | None:
    """Compatibility entry point; selected scaffolds must provide explicit facts."""
    if facts is None and required:
        raise ValueError("numbering_scaffold facts required for selected scaffold")
    parent = {"numbering_scaffold": facts} if facts is not None else {}
    return plan_from_parent(chain, parent)


def effective_sub_locant(plan: NumberingPlan, atom: int) -> int | None:
    """Plain int matching legacy _sub_locant (bridgehead → chain index + 1)."""
    lab = locant(plan, atom)
    if lab is None:
        return None
    if lab.isdigit():
        return int(lab)
    try:
        return plan.atom_order.index(atom) + 1
    except ValueError:
        return None
