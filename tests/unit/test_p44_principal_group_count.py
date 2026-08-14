# IUPAC: P-44.1.1
# Layer: L2
"""Principal-group count precedes later parent-selection criteria."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
import hashlib
import namepredict.namer as namer_module
from namepredict.namer import SMILESNNamer
from namepredict.layer2 import kind_registry
from namepredict.layer2.parent_candidate import from_parent_dict, principal_contract_kind
from namepredict.layer2.scoring import _score_parent
from namepredict.types import NameResult

COUNT_CASES = [
    ({"kind": "alcohol", "principal_group_count": 2}, (5, 2)),
    ({"kind": "amine", "principal_group_count": 3}, (3, 3)),
    ({"kind": "ketone", "principal_group_count": 2}, (6, 2)),
]


@pytest.mark.parametrize("parent,expected", COUNT_CASES)
def test_typed_principal_group_contract(parent, expected):
    facts = from_parent_dict(parent).facts
    assert (facts.principal_group_class, facts.principal_group_count) == expected


def test_scoring_starts_with_typed_p44_facts():
    parent = {"kind": "alcohol", "principal_group_count": 2}
    assert _score_parent({}, parent)[:2] == (5, 2)


def test_higher_class_precedes_larger_lower_class_count():
    acid = {"kind": "acid", "principal_group_count": 1}
    triol = {"kind": "alcohol", "principal_group_count": 3}
    amine = {"kind": "amine", "principal_group_count": 1}
    assert _score_parent({}, acid) > _score_parent({}, triol) > _score_parent({}, amine)


def test_missing_principal_count_is_not_silently_one():
    with pytest.raises(ValueError, match="principal_group_count"):
        from_parent_dict({"kind": "acid"})


def _result(en):
    return NameResult(en=en, zh=en, success=True, source="iupac", time_ms=0.0)


def _prepared(kind, complete):
    return ({"kind": kind}, [], complete)


def test_higher_phase_partial_preserves_seniority(monkeypatch):
    phases = [[{"kind": "acid"}], [{"kind": "alcohol"}]]
    monkeypatch.setattr(namer_module, "_candidate_phases", lambda info, depth: phases)
    monkeypatch.setattr(
        namer_module, "_prepare_candidate",
        lambda info, parent, **kwargs: _prepared(parent["kind"], parent["kind"] == "alcohol"),
    )
    monkeypatch.setattr(namer_module, "_assemble_candidate", lambda p, *a, **k: _result(p["kind"]))
    assert namer_module._run_candidates({}, depth=0, t0=0).en == "acid"


def test_higher_phase_complete_never_downgrades(monkeypatch):
    phases = [[{"kind": "acid"}], [{"kind": "alcohol"}]]
    monkeypatch.setattr(namer_module, "_candidate_phases", lambda info, depth: phases)
    monkeypatch.setattr(
        namer_module, "_prepare_candidate",
        lambda info, parent, **kwargs: _prepared(parent["kind"], True),
    )
    monkeypatch.setattr(namer_module, "_assemble_candidate", lambda p, *a, **k: _result(p["kind"]))
    assert namer_module._run_candidates({}, depth=0, t0=0).en == "acid"

CASES = [
    ("OCCC(CCCl)C(O)C", "3-(2-chloroethyl)pentane-1,4-diol", None),
    ("OCCC(CCCCl)C(O)C", "3-(3-chloropropyl)pentane-1,4-diol", None),
    ("OCCC(CCCCCl)C(O)C", "3-(4-chlorobutyl)pentane-1,4-diol", None),
    ("OCCC(CCCCBr)C(O)C", "3-(4-bromobutyl)pentane-1,4-diol", None),
    # Near-neighbour negative: the already claimable C1 arm stays unchanged.
    ("OCCC(CCl)C(O)C", "3-(chloromethyl)pentane-1,4-diol", None),
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_principal_group_count_controls_entry_parent(smiles, en, zh):
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(result.zh) == normalize_zh(zh)
