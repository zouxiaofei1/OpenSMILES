# IUPAC: P-14
# Layer: L4
"""NumberingPlan pure data + locant / locant_int API (no engine)."""
from __future__ import annotations

from namepredict.layer4.locants import NumberingPlan, locant, locant_int, make_plan


def _simple_plan() -> NumberingPlan:
    return make_plan(
        scaffold_id="test",
        atom_order=[0, 1, 2],
        labels=["1", "2", "3"],
        sub_atoms=frozenset({1, 2}),
        constraints_applied=("principal",),
    )


def _bridge_plan() -> NumberingPlan:
    return make_plan(
        scaffold_id="fused",
        atom_order=[10, 11, 12, 13, 14],
        labels=["1", "2", "3", "3a", "4"],
        sub_atoms=frozenset({10, 11, 12, 14}),
        constraints_applied=(),
    )


def test_locant_string_lookup():
    plan = _simple_plan()
    assert locant(plan, 1) == "2"
    assert locant(plan, 0) == "1"
    assert locant(plan, 2) == "3"


def test_locant_int_plain():
    """Plain labels scale as n*10 (e.g. '2' → 20)."""
    plan = _simple_plan()
    assert locant_int(plan, 1) == 20
    assert locant_int(plan, 0) == 10


def test_locant_unknown_atom_none():
    plan = _simple_plan()
    assert locant(plan, 99) is None
    assert locant_int(plan, 99) is None


def test_label_to_atom_roundtrip():
    plan = _simple_plan()
    assert plan.label_to_atom["2"] == 1
    assert plan.atom_to_label[1] == "2"
    for atom, label in plan.atom_to_label.items():
        assert plan.label_to_atom[label] == atom


def test_bridgehead_3a_total_order():
    """Scaled-int convention: '3'→30, '3a'→31, '4'→40 so 3 < 3a < 4."""
    plan = _bridge_plan()
    assert locant(plan, 13) == "3a"
    i3 = locant_int(plan, 12)
    i3a = locant_int(plan, 13)
    i4 = locant_int(plan, 14)
    assert i3 == 30 and i3a == 31 and i4 == 40
    assert i3 < i3a < i4


def test_plan_fields():
    plan = _simple_plan()
    assert plan.scaffold_id == "test"
    assert plan.atom_order == (0, 1, 2)
    assert plan.labels == ("1", "2", "3")
    assert plan.sub_atoms == frozenset({1, 2})
    assert plan.constraints_applied == ("principal",)
