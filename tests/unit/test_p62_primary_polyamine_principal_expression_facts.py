# IUPAC: P-62.2.1
# Layer: L1,L2,L3,L4,L5
import pytest
from rdkit import Chem

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.principal_expression import PrincipalChargeState, PrincipalRelation
from namepredict.namer import SMILESNNamer


CASES = [
    ("CCN", "amine", 1, "ethanamine", "乙胺"),
    ("CC(N)C", "amine", 1, "propan-2-amine", "丙-2-胺"),
    ("NCCN", "amine", 2, "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("NCCCN", "amine", 2, "propane-1,3-diamine", "丙烷-1,3-二胺"),
    ("NCC(N)CN", "amine", 3, "propane-1,2,3-triamine", "丙烷-1,2,3-三胺"),
]


def _typed_parent(smiles, kind, count):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("kind") == kind
                and p.get("principal_expression_facts")
                and p["principal_expression_facts"].multiplicity == count)


@pytest.mark.parametrize("smiles,kind,count,en,zh", CASES)
def test_primary_polyamine_typed_facts_and_names(smiles, kind, count, en, zh):
    facts = _typed_parent(smiles, kind, count)["principal_expression_facts"]
    assert facts.group_class.value == "amine"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert len(facts.attachment_atoms) == count
    expected = {a["c_idx"] for a in analyze(Chem.MolFromSmiles(smiles))["amines"]}
    assert set(facts.attachment_atoms) == expected
    assert facts.charge_state is PrincipalChargeState.NEUTRAL
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["CCNCC", "CN(C)C"])
def test_sec_tert_amine_single_typed_path(smiles):
    """二级/三级胺走单个胺 typed 路径（母体含 N 链 + N- 取代基），不产生歧义多候选。"""
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    amines = [p for p in parents if p.get("kind") == "amine" and p.get("principal_expression_facts")]
    assert amines
    assert all(p["principal_expression_facts"].multiplicity == 1 for p in amines)
    result = SMILESNNamer().name(smiles)
    assert result.success


def test_aniline_stays_on_retained_boundary():
    for smiles, en in [("Nc1ccccc1", "aniline"), ("COc1ccc(N)cc1", "4-methoxyaniline")]:
        parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
        assert any(p.get("kind") == "amine" and p.get("scaffold_id") == "benzene" for p in parents)
        result = SMILESNNamer().name(smiles)
        assert normalize_en(result.en) == en


def test_amino_acid_keeps_acid_principal():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("NCC(=O)O")))
    facts = [p.get("principal_expression_facts") for p in parents]
    assert not any(f and f.group_class.value == "amine" for f in facts)
