# IUPAC: P-63.1.1 / P-63.1.2
# Layer: L2,L4,L5
import pytest
from rdkit import Chem

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.principal_expression import PrincipalChargeState, PrincipalRelation
from namepredict.namer import SMILESNNamer


CASES = [
    ("CCCO", "alcohol", 1, "propan-1-ol", "丙-1-醇"),
    ("CC(O)C", "alcohol", 1, "propan-2-ol", "丙-2-醇"),
    ("OCCCO", "alcohol", 2, "propane-1,3-diol", "丙烷-1,3-二醇"),
    ("CC(O)C(O)C", "alcohol", 2, "butane-2,3-diol", "丁烷-2,3-二醇"),
    ("OCC(O)CO", "alcohol", 3, "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("OC1CCCCC1", "alcohol", 1, "cyclohexanol", "环己醇"),
    ("OC1CCC(O)CC1", "alcohol", 2, "cyclohexane-1,4-diol", "环己烷-1,4-二醇"),
]


def _typed_parent(smiles: str, kind: str, count: int):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("kind") == kind
                and p.get("principal_expression_facts")
                and p["principal_expression_facts"].multiplicity == count)


@pytest.mark.parametrize("smiles,kind,count,en,zh", CASES)
def test_alcohol_principal_expression_facts_and_names(smiles, kind, count, en, zh):
    parent = _typed_parent(smiles, kind, count)
    facts = parent["principal_expression_facts"]
    assert facts.group_class.value == "alcohol"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert len(facts.attachment_atoms) == count
    assert facts.charge_state is PrincipalChargeState.NEUTRAL
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_hydroxy_acid_keeps_acid_as_principal_group():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CC(O)C(=O)O")))
    facts = [p.get("principal_expression_facts") for p in parents]
    assert not any(f and f.group_class.value == "alcohol" for f in facts)
    result = SMILESNNamer().name("CC(O)C(=O)O")
    assert result.success
    assert normalize_en(result.en) == "2-hydroxypropanoic acid"
    assert normalize_zh(result.zh) == "2-羟基丙酸"


def test_ether_has_no_alcohol_principal_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCOCC")))
    facts = [p.get("principal_expression_facts") for p in parents]
    assert not any(f and f.group_class.value == "alcohol" for f in facts)
