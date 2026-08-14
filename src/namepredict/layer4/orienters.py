"""Orientation helpers still used by layer4 locant computation.

The kind-dispatched orienters (`_kind_orienters` and friends) were replaced by
the candidate-based numbering engine (`numbering_engine.orient_numbering`);
only `_typed_group_atoms` survives for FG locant computation in locant_calc.
"""
from __future__ import annotations


def _typed_group_atoms(parent: dict, group: str) -> list[int]:
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts and facts.group_class.value == group else []
