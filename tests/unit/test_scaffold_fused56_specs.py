# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4
"""Fused56 ScaffoldSpec registry + labels authority + thin parent mapping."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer2.specs import (
    FUSED56_LABELS,
    FUSED56_SPECS,
    ScaffoldSpec,
    all_specs,
    get_spec,
    numbering_scaffold_facts,
)
from namepredict.layer4.locants.adapt import plan_from_chain
from namepredict.namer import SMILESNNamer

# Labels authority: 1..3a..7a (9 roles).
_EXPECTED_LABELS = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")

# Positive: fused56 retained parents must stay registered + name-stable.

# Negative: non-fused56 parents must not be fused56 specs.
NEG_KINDS = ("benzene", "naphthalene", "cycloalkane", "furan")


@pytest.mark.parametrize(
    "sid",
    [
        "benzofuran", "benzofuranamine",
        "benzothiophene", "benzothiophenol",
        "benzothiazole", "benzothiazolamine",
        "benzoxazole", "benzoxazolamine",
    ],
)
def test_fused56_specs_registered(sid: str) -> None:
    sp = get_spec(sid)
    assert isinstance(sp, ScaffoldSpec)
    assert sp.id == sid
    assert sp.naming_class == "fused56"
    assert sp.n_rings == 2
    assert sp.ring == "hetero"
    assert sp.retained is True
    assert sp.numbering.standard_path == _EXPECTED_LABELS
    assert sp.numbering.mode in ("fused56_fixed", "fixed_roles")
    assert sp.id in {s.id for s in FUSED56_SPECS}
    assert get_spec(sid) is not None


def test_fused56_labels_constant() -> None:
    assert FUSED56_LABELS == _EXPECTED_LABELS
    assert len(FUSED56_LABELS) == 9
    assert "3a" in FUSED56_LABELS and "7a" in FUSED56_LABELS


def test_all_specs_includes_fused56() -> None:
    ids = {s.id for s in all_specs()}
    assert "benzofuran" in ids
    assert "benzothiazole" in ids
    assert "cycloalkane" in ids  # carbocycle preserved


@pytest.mark.parametrize("kind", NEG_KINDS)
def test_non_fused56_not_wrong_parent_spec(kind: str) -> None:
    sp = get_spec(kind)
    if sp is None:
        return
    assert sp.naming_class != "fused56"


def test_plan_from_chain_uses_spec_labels() -> None:
    chain = list(range(9))
    plan = plan_from_chain(chain, "benzofuran", numbering_scaffold_facts("benzofuran", len(chain)))
    assert plan is not None
    assert plan.labels == _EXPECTED_LABELS
    assert plan.scaffold_id == "benzofuran"


def test_plan_from_chain_benzene_none() -> None:
    assert plan_from_chain(list(range(6)), "benzene", numbering_scaffold_facts("benzene", 6)) is None


def test_l2_fused56_labels_are_parent_plan_authority() -> None:
    facts = numbering_scaffold_facts("benzofuran", len(_EXPECTED_LABELS))
    assert facts is not None
    assert FUSED56_LABELS == facts["labels"] == _EXPECTED_LABELS
