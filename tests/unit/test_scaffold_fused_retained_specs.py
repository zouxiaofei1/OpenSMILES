# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4
"""Retained fused aza scaffolds (indole/bim/quinoline/naph) as ScaffoldSpec."""
from __future__ import annotations

import pytest

from namepredict.layer2.scaffold.specs import (
    FUSED56_LABELS,
    NAPH_LABELS,
    ScaffoldSpec,
    fused56_kind_ids,
    get_spec,
    kind_ids_for,
    naph_kind_ids,
    numbering_scaffold_facts,
)
from namepredict.layer4.locants.adapt import plan_from_chain

_F56_LABELS = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")
_NAPH_LABELS = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")

# Negative: wrong naming_class / not these families
NEG_CASES = [
    ("benzene", None),  # no retained fused Spec
    ("furan", None),
    ("cycloalkane", "carbocycle"),  # carbo free, not fused56/naph
]


@pytest.mark.parametrize(
    "sid,nclass,n_labs",
    [
        ("indole", "fused56", 9),
        ("indazole", "fused56", 9),
        ("benzimidazole", "fused56", 9),
        ("benzimidazolamine", "fused56", 9),
        ("quinoline", "naph_family", 10),
        ("isoquinoline", "naph_family", 10),
        ("naphthalene", "naph_family", 10),
        ("naphthalenol", "naph_family", 10),
        ("naphthalenediol", "naph_family", 10),
    ],
)
def test_retained_fused_specs_registered(sid: str, nclass: str, n_labs: int) -> None:
    sp = get_spec(sid)
    assert isinstance(sp, ScaffoldSpec)
    assert sp.id == sid
    assert sp.naming_class == nclass
    assert sp.n_rings == 2
    assert sp.retained is True
    assert len(sp.numbering.standard_path) == n_labs
    if nclass == "fused56":
        assert sp.numbering.mode == "fused56_fixed"
        assert sp.numbering.standard_path == _F56_LABELS
        assert sp.ring == "hetero"
    else:
        assert sp.numbering.mode in ("naph_family", "naph_fixed")
        assert sp.numbering.standard_path == _NAPH_LABELS


def test_labels_constants() -> None:
    assert FUSED56_LABELS == _F56_LABELS
    assert NAPH_LABELS == _NAPH_LABELS
    assert len(NAPH_LABELS) == 10
    assert "4a" in NAPH_LABELS and "8a" in NAPH_LABELS


def test_kind_id_helpers() -> None:
    f56 = fused56_kind_ids()
    assert "indole" in f56 and "benzimidazole" in f56
    assert "benzofuran" in f56  # prior fused56 retained
    assert "quinoline" not in f56
    naph = naph_kind_ids()
    assert "quinoline" in naph and "naphthalene" in naph
    assert "indole" not in naph
    assert kind_ids_for("fused56") == f56
    assert kind_ids_for("naph_family") == naph


@pytest.mark.parametrize("kind,nclass", NEG_CASES)
def test_negative_not_wrong_class(kind: str, nclass: str | None) -> None:
    sp = get_spec(kind)
    if nclass is None:
        assert sp is None or sp.naming_class not in ("fused56", "naph_family")
        return
    assert sp is not None
    assert sp.naming_class == nclass
    assert sp.naming_class not in ("fused56", "naph_family") or kind == "never"


def test_plan_from_chain_indole_9() -> None:
    plan = plan_from_chain(list(range(9)), "indole", numbering_scaffold_facts("indole", 9))
    assert plan is not None
    assert plan.labels == _F56_LABELS
    assert plan.scaffold_id == "indole"


def test_plan_from_chain_quinoline_10() -> None:
    plan = plan_from_chain(list(range(10)), "quinoline", numbering_scaffold_facts("quinoline", 10))
    assert plan is not None
    assert plan.labels == _NAPH_LABELS
    assert plan.scaffold_id == "quinoline"


def test_plan_from_chain_naphthalene_10() -> None:
    plan = plan_from_chain(list(range(10)), "naphthalene", numbering_scaffold_facts("naphthalene", 10))
    assert plan is not None
    assert plan.labels == _NAPH_LABELS
