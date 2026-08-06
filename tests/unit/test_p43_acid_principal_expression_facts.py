# IUPAC: P-43.1
# Layer: L2,L4,L5
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.principal_expression import PrincipalChargeState, PrincipalRelation, express_chain_principal
from namepredict.layer2.principal_parent import select_principal_parent_skeletons
from namepredict.namer import SMILESNNamer


CASES = [
    ("CCCC(=O)O", "acid", 1, PrincipalRelation.IN_SKELETON, "neutral", "butanoic acid", "丁酸"),
    ("O=C(O)CCC(=O)O", "diacid", 2, PrincipalRelation.IN_SKELETON, "neutral", "butanedioic acid", "丁二酸"),
    ("O=C(O)C1CCCCC1", "cycloalkane", 1, PrincipalRelation.EXOCYCLIC, "neutral", "cyclohexanecarboxylic acid", "环己烷甲酸"),
    ("CCCC(=O)[O-]", "acid", 1, PrincipalRelation.IN_SKELETON, "anion", "butanoate", "丁酸根"),
]


def _parent(smiles, kind):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(parent for parent in parents if parent.get("kind") == kind)


@pytest.mark.parametrize("smiles,kind,count,relation,charge,en,zh", CASES)
def test_acid_principal_expression_facts(smiles, kind, count, relation, charge, en, zh):
    facts = _parent(smiles, kind)["principal_expression_facts"]
    assert facts.group_class.value == "acid"
    assert facts.multiplicity == count
    assert facts.relation is relation
    assert len(facts.occurrence_ids) == count
    assert facts.characteristic_atoms
    assert facts.attachment_atoms
    assert facts.charge_state == charge
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_ester_near_neighbor_has_no_acid_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCCC(=O)OC")))
    ester = next(parent for parent in parents if parent.get("kind") == "ester")
    assert "principal_expression_facts" not in ester


def test_partially_deprotonated_diacid_has_mixed_charge_state():
    parent = _parent("O=C([O-])CCC(=O)O", "diacid")
    assert parent["principal_expression_facts"].charge_state is PrincipalChargeState.MIXED


def test_facts_include_only_occurrences_covered_by_selected_skeleton():
    info = analyze(Chem.MolFromSmiles("C(C(=O)O)(C(=O)O)C(=O)O"))
    selected = select_principal_parent_skeletons(info)
    skeleton = selected.skeletons.candidates[0]
    parent = express_chain_principal(info, selected.principal, skeleton)
    facts = parent["principal_expression_facts"]
    assert facts.multiplicity == len(skeleton.covered_principal_ids) == 2
    assert set(facts.occurrence_ids) == set(skeleton.covered_principal_ids)
