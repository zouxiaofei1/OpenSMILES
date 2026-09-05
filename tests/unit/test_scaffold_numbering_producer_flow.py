# IUPAC: P-14 / P-25
# Layer: L2,L4
"""Selected retained parents must carry facts before L4 numbering."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.ring_scaffold import all_specs
from namepredict.layer4.numbering import number
from namepredict.namer import SMILESNNamer


_CASES = [
    ("c1ccc2[nH]ccc2c1", "indole", "1H-indole"),
    ("c1ccc2ccccc2c1", "naphthalene", "naphthalene"),
    ("C1CCCCC1", "cycloalkane", "cyclohexane"),
]


def _selected(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return select_parent(analyze(mol))


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
