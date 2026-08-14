# IUPAC: P-14 / P-25
# Layer: L4
"""NumberingPlan / plan_from_chain contract tests.

正交化后 _ALL_SPECS 的 standard_path 全为空，保留 scaffold 的固定标签（3a/7a/…）
已废弃，plan 机制仅在有显式 facts 时产出。这里测纯 plan 机制契约，不依赖具体 scaffold。
"""
from __future__ import annotations

import pytest

from namepredict.layer4.locants.adapt import (
    effective_sub_locant,
    plan_from_chain,
    plan_from_parent,
)
from namepredict.layer4.locants.plan import locant, make_plan


def _labels(atom_order: tuple, labs: tuple):
    return make_plan("fused", list(atom_order), labs)


def test_make_plan_maps_labels_both_ways() -> None:
    plan = _labels((0, 1, 2), ("1", "3a", "7a"))
    assert plan.atom_order == (0, 1, 2)
    assert plan.atom_to_label == {0: "1", 1: "3a", 2: "7a"}
    assert plan.label_to_atom == {"1": 0, "3a": 1, "7a": 2}


def test_locant_queries_atom_label() -> None:
    plan = _labels((0, 1, 2), ("1", "3a", "7a"))
    assert locant(plan, 1) == "3a"
    assert locant(plan, 9) is None


def test_effective_sub_locant_numeric_label_returns_int() -> None:
    plan = _labels((0, 1, 2), ("1", "2", "3"))
    assert effective_sub_locant(plan, 1) == 2


def test_effective_sub_locant_fused_label_uses_index_plus_one() -> None:
    # 非数字 label（3a/7a）→ 回落 chain.index + 1（bridgehead 行为）。
    plan = _labels((100, 101, 102), ("1", "3a", "7a"))
    assert effective_sub_locant(plan, 101) == 2


def test_effective_unknown_atom_none() -> None:
    plan = _labels((0, 1, 2), ("1", "3a", "7a"))
    assert effective_sub_locant(plan, 999) is None


def test_plan_from_chain_none_without_facts() -> None:
    assert plan_from_chain([0, 1, 2], "alkane") is None
    assert plan_from_chain([0, 1, 2], "naphthalene") is None


def test_plan_from_chain_required_without_facts_raises() -> None:
    with pytest.raises(ValueError):
        plan_from_chain([0, 1, 2], "naphthalene", required=True)


def test_plan_from_chain_with_facts_builds_plan() -> None:
    facts = {"scaffold_id": "fused", "labels": ("1", "3a", "7a")}
    plan = plan_from_chain([0, 1, 2], "fused", facts)
    assert plan is not None
    assert plan.labels == ("1", "3a", "7a")
    assert plan.scaffold_id == "fused"


def test_plan_from_parent_wrong_chain_length_none() -> None:
    facts = {"scaffold_id": "fused", "labels": ("1", "2", "3")}
    assert plan_from_parent([0, 1], {"numbering_scaffold": facts}) is None


def test_plan_from_parent_non_mapping_facts_none() -> None:
    assert plan_from_parent([0, 1, 2], {"numbering_scaffold": None}) is None
    assert plan_from_parent([0, 1, 2], {"numbering_scaffold": object()}) is None


def test_plan_from_chain_non_numeric_labels_rejected() -> None:
    # labels 必须全部为 str；非 str 拒绝 → None。
    facts = {"scaffold_id": "fused", "labels": (1, 2, 3)}
    assert plan_from_chain([0, 1, 2], "fused", facts) is None
