# IUPAC: P-14
# Layer: L4
"""NumberingPlan pure data + locant API."""
from __future__ import annotations

from namepredict.layer4.locants import NumberingPlan, locant, make_plan


def _simple_plan() -> NumberingPlan:
    return make_plan(
        scaffold_id="test",
        atom_order=[0, 1, 2],
        labels=["1", "2", "3"],
        sub_atoms=frozenset({1, 2}),
        constraints_applied=("principal",),
    )


def test_locant_string_lookup():
    plan = _simple_plan()
    assert locant(plan, 1) == "2"
    assert locant(plan, 0) == "1"
    assert locant(plan, 2) == "3"


def test_locant_unknown_atom_none():
    plan = _simple_plan()
    assert locant(plan, 99) is None


def test_label_to_atom_roundtrip():
    plan = _simple_plan()
    assert plan.label_to_atom["2"] == 1
    assert plan.atom_to_label[1] == "2"
    for atom, label in plan.atom_to_label.items():
        assert plan.label_to_atom[label] == atom


def test_plan_fields():
    plan = _simple_plan()
    assert plan.scaffold_id == "test"
    assert plan.atom_order == (0, 1, 2)
    assert plan.labels == ("1", "2", "3")
    assert plan.sub_atoms == frozenset({1, 2})
    assert plan.constraints_applied == ("principal",)
