# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4
"""Fused56 ScaffoldSpec registry + labels authority + thin parent mapping."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.benzofuran import _try_benzofuran_parent
from namepredict.layer2.scaffold.benzothiazole import _try_benzothiazole_parent
from namepredict.layer2.scaffold.builders.fused56 import (
    FUSED56_LABELS,
    scaffold_id_for_kind,
)
from namepredict.layer2.scaffold.specs import (
    FUSED56_SPECS,
    ScaffoldSpec,
    all_specs,
    get_spec,
)
from namepredict.layer2.scaffold.specs import numbering_scaffold_facts
from namepredict.layer4.locants.adapt import plan_from_chain
from namepredict.namer import SMILESNNamer

# Labels authority: 1..3a..7a (9 roles).
_EXPECTED_LABELS = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")

# Positive: fused56 retained parents must stay registered + name-stable.
POS_CASES = [
    ("c1ccc2occc2c1", "benzofuran", "benzofuran", "苯并呋喃"),
    ("c1ccc2sccc2c1", "benzothiophene", "1-benzothiophene", "苯并[b]噻吩"),
    ("c1ccc2scnc2c1", "benzothiazole", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("c1ccc2ocnc2c1", "benzoxazole", "1,3-benzoxazole", "1,3-苯并噁唑"),
    ("O1C(=CC2=C1C=CC=C2)N", "benzofuranamine", "benzofuran-2-amine", "苯并呋喃-2-胺"),
]

# Negative: non-fused56 parents must not be fused56 specs.
NEG_KINDS = ("benzene", "naphthalene", "cycloalkane", "furan")


def _mol(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    return mol


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


@pytest.mark.parametrize("smiles,kind,en,zh", POS_CASES)
def test_e2e_fused56_names_stable(smiles, kind, en, zh) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_parent_carries_scaffold_id() -> None:
    parent = _try_benzofuran_parent(analyze(_mol("c1ccc2occc2c1")))
    assert parent is not None
    assert parent.get("kind") == "benzofuran"
    assert parent.get("scaffold_id") == "benzofuran"
    assert len(parent.get("chain") or []) == 9


def test_btz_parent_scaffold_id() -> None:
    parent = _try_benzothiazole_parent(analyze(_mol("c1ccc2scnc2c1")))
    assert parent is not None
    assert parent.get("scaffold_id") == "benzothiazole"


def test_plan_from_chain_uses_spec_labels() -> None:
    chain = list(range(9))
    plan = plan_from_chain(chain, "benzofuran", numbering_scaffold_facts("benzofuran", len(chain)))
    assert plan is not None
    assert plan.labels == _EXPECTED_LABELS
    assert plan.scaffold_id == "benzofuran"


def test_plan_from_chain_benzene_none() -> None:
    assert plan_from_chain(list(range(6)), "benzene", numbering_scaffold_facts("benzene", 6)) is None


def test_scaffold_id_helpers() -> None:
    assert scaffold_id_for_kind("benzofuran") == "benzofuran"
    assert scaffold_id_for_kind("benzothiazolamine") == "benzothiazolamine"
    assert scaffold_id_for_kind("benzene") is None


def test_l2_fused56_labels_are_parent_plan_authority() -> None:
    facts = numbering_scaffold_facts("benzofuran", len(_EXPECTED_LABELS))
    assert facts is not None
    assert FUSED56_LABELS == facts["labels"] == _EXPECTED_LABELS
