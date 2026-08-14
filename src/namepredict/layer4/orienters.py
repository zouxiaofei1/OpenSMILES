from __future__ import annotations


def _typed_group_atoms(parent: dict, group: str) -> list[int]:
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts and facts.group_class.value == group else []
