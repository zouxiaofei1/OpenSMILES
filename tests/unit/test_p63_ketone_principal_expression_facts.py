# IUPAC: P-63.2.1
# Layer: L2,L4,L5
import pytest
from rdkit import Chem

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.principal_expression import PrincipalChargeState, PrincipalRelation
from namepredict.namer import SMILESNNamer


CASES = [
    ("CC(=O)CC", "ketone", 1, "butan-2-one", "丁-2-酮"),
    ("CCC(=O)CC", "ketone", 1, "pentan-3-one", "戊-3-酮"),
    ("CC(=O)CC(=O)C", "ketone", 2, "pentane-2,4-dione", "戊-2,4-二酮"),
    ("CC(=O)C(C)=O", "ketone", 2, "butane-2,3-dione", "丁-2,3-二酮"),
]


def _parent(smiles, kind):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(parent for parent in parents if parent.get("kind") == kind)


@pytest.mark.parametrize("smiles,kind,count,en,zh", CASES)
def test_ketone_principal_expression_facts_and_names(smiles, kind, count, en, zh):
    facts = _parent(smiles, kind)["principal_expression_facts"]
    assert facts.group_class.value == "ketone"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert facts.characteristic_atoms and facts.attachment_atoms
    assert facts.charge_state is PrincipalChargeState.NEUTRAL
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_alcohol_near_neighbor_has_no_ketone_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCCCO")))
    alcohol = next(parent for parent in parents if parent.get("kind") == "alcohol")
    facts = alcohol.get("principal_expression_facts")
    assert facts is None or facts.group_class.value != "ketone"
