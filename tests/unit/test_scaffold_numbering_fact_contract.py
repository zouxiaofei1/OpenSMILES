# IUPAC: P-14 / P-25
# Layer: L2/L4 contract
"""Parent-owned numbering facts are the sole retained-scaffold label authority."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.layer2.specs import all_specs, numbering_scaffold_facts
from namepredict.layer4.locants.adapt import plan_from_parent


_SUPPORTED = tuple(
    sp for sp in all_specs()
    if sp.numbering.standard_path and sp.numbering.materialize_plan
)


@pytest.mark.parametrize("spec", _SUPPORTED, ids=lambda sp: sp.id)
def test_parent_numbering_facts_produce_the_spec_labels(spec) -> None:
    chain = list(range(len(spec.numbering.standard_path)))
    facts = numbering_scaffold_facts(spec.id, len(chain))
    assert facts == {
        "scaffold_id": spec.id,
        "labels": spec.numbering.standard_path,
        "relative_stereo": None,
    }
    plan = plan_from_parent(chain, {"numbering_scaffold": facts})
    assert plan is not None
    assert plan.scaffold_id == spec.id
    assert plan.atom_order == tuple(chain)
    assert plan.labels == spec.numbering.standard_path


@pytest.mark.parametrize("size", range(3, 11))
def test_cyclo_dynamic_parent_facts_produce_numeric_labels(size: int) -> None:
    chain = list(range(size))
    facts = numbering_scaffold_facts("cycloalkane", size)
    assert facts is not None
    assert facts["labels"] == tuple(str(i) for i in range(1, size + 1))
    plan = plan_from_parent(chain, {"numbering_scaffold": facts})
    assert plan is not None
    assert plan.labels == facts["labels"]


@pytest.mark.parametrize("spec_id,size", [
    ("benzene", 6), ("quinazoline", 10), ("unknown", 9),
    ("cycloalkane", 2), ("cycloalkane", 11),
])
def test_unsupported_or_unknown_scaffold_has_no_numbering_plan(spec_id: str, size: int) -> None:
    facts = numbering_scaffold_facts(spec_id, size)
    assert facts is None
    assert plan_from_parent(list(range(size)), {"numbering_scaffold": facts}) is None


def test_l4_adapter_does_not_import_layer2_or_keep_kind_label_whitelists() -> None:
    path = Path(__file__).parents[2] / "src/namepredict/layer4/locants/adapt.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imports = [node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert not any(name.startswith("namepredict.layer2") for name in imports)
    source = path.read_text(encoding="utf-8")
    assert "FUSED56_KINDS" not in source
    assert "Q_KINDS" not in source
    assert "NAPH_KINDS" not in source
    assert "ANTHRA_KINDS" not in source
