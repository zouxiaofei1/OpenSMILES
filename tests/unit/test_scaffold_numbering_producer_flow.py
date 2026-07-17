# IUPAC: P-14 / P-25
# Layer: L2,L4
"""Selected retained parents must carry facts before L4 numbering."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.scaffold.specs import all_specs
from namepredict.layer4.numbering import number
from namepredict.namer import SMILESNNamer


_CASES = [
    ("c1ccc2[nH]ccc2c1", "indole", "1H-indole"),
    ("c1ccc2ncccc2c1", "quinoline", "quinoline"),
    ("c1ccc2ccccc2c1", "naphthalene", "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1", "anthracene", "anthracene"),
    ("C1CCCCC1", "cycloalkane", "cyclohexane"),
]


def _selected(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return select_parent(analyze(mol))


@pytest.mark.parametrize("smiles,kind,_", _CASES)
def test_selected_parent_materializes_facts_that_l4_uses(smiles: str, kind: str, _: str) -> None:
    parent = _selected(smiles)
    assert parent["kind"] == kind
    facts = parent.get("numbering_scaffold")
    assert facts is not None
    numbered = number(parent, [])
    assert numbered["parent"]["numbering"].labels == facts["labels"]


def test_all_plan_enabled_specs_have_a_label_contract() -> None:
    enabled = [sp for sp in all_specs() if sp.numbering.materialize_plan and sp.numbering.standard_path]
    assert {sp.id for sp in enabled} == {
        "anthracene", "anthraquinone", "benzofuran", "benzofuranamine",
        "benzothiazole", "benzothiazolamine", "benzothiophene", "benzothiophenol",
        "benzimidazole", "benzimidazolamine", "benzoxazole", "benzoxazolamine",
        "indazole", "indazolecarbaldehyde", "indazolecarbonitrile", "indole",
        "indolecarboxylic", "isoquinoline", "naphthalene", "naphthalenecarboxylic",
        "quinoline", "quinolinecarboxylic", "quinolinol",
    }


def test_scaffold_parent_without_materialized_facts_is_rejected() -> None:
    parent = _selected("c1ccc2[nH]ccc2c1")
    parent.pop("numbering_scaffold")
    with pytest.raises(ValueError, match="numbering_scaffold"):
        number(parent, [])


@pytest.mark.parametrize("smiles,_,expected", _CASES)
def test_representative_public_names_are_unchanged(smiles: str, _: str, expected: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert result.en == expected
