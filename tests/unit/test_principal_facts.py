# 合并自 5 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_p43_acid_principal_expression_facts.py: 
test_p62_neutral_amine_arm_facts.py: 
test_p62_primary_polyamine_principal_expression_facts.py: 
test_p63_alcohol_principal_expression_facts.py: 
test_p63_ketone_principal_expression_facts.py: 
"""
from __future__ import annotations

import pytest

from opensmiles.layer1.analyzer import analyze
from opensmiles.layer1.functional_group_inventory import FunctionalGroupClass as FG, inventory_from_info
from opensmiles.layer2.parent_select import _collect_candidates
from opensmiles.layer2.principal_expression import PrincipalRelation, express_chain_principal
from opensmiles.layer2.parent_select import select_principal_parent_skeletons
from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# 合并自 test_p43_acid_principal_expression_facts.py
# IUPAC: P-43.1
# Layer: L2,L4,L5
# ==========================================================================
p43_acid_principal_expression_facts__CASES = [
    ("CCCC(=O)O", "acid", 1, PrincipalRelation.IN_SKELETON, "butanoic acid", "丁酸"),
    ("O=C(O)CCC(=O)O", "acid", 2, PrincipalRelation.IN_SKELETON, "butanedioic acid", "丁二酸"),
    ("CCCC(=O)[O-]", "acid", 1, PrincipalRelation.IN_SKELETON, "butanoate", "丁酸根"),
]


def p43_acid_principal_expression_facts___parent(smiles, kind, count):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(parent for parent in parents
                if parent.get("kind") == kind
                and parent.get("principal_expression_facts")
                and parent["principal_expression_facts"].multiplicity == count)


@pytest.mark.parametrize("smiles,kind,count,relation,en,zh", p43_acid_principal_expression_facts__CASES)
def test_acid_principal_expression_facts(smiles, kind, count, relation, en, zh):
    facts = p43_acid_principal_expression_facts___parent(smiles, kind, count)["principal_expression_facts"]
    assert facts.group_class.value == "acid"
    assert facts.multiplicity == count
    assert facts.relation is relation
    assert len(facts.occurrence_ids) == count
    assert facts.characteristic_atoms
    assert facts.attachment_atoms
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_facts_include_only_occurrences_covered_by_selected_skeleton():
    info = analyze(Chem.MolFromSmiles("C(C(=O)O)(C(=O)O)C(=O)O"))
    selected = select_principal_parent_skeletons(info)
    skeleton = selected.skeletons.candidates[0]
    parent = express_chain_principal(info, selected.principal, skeleton)
    facts = parent["principal_expression_facts"]
    assert facts.multiplicity == len(skeleton.covered_principal_ids) == 2
    assert set(facts.occurrence_ids) == set(skeleton.covered_principal_ids)


# ==========================================================================
# 合并自 test_p62_neutral_amine_arm_facts.py
# IUPAC: P-62.2.2.1
# Layer: L2,L3,L4,L5
# ==========================================================================
def p62_neutral_amine_arm_facts___arm_parent(smiles):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("neutral_amine_arm_facts"))


def test_primary_amine_has_no_n_arm_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCN")))
    assert not any(p.get("neutral_amine_arm_facts") for p in parents)


# ==========================================================================
# 合并自 test_p62_primary_polyamine_principal_expression_facts.py
# IUPAC: P-62.2.1
# Layer: L1,L2,L3,L4,L5
# ==========================================================================
p62_primary_polyamine_principal_expression_facts__CASES = [
    ("CCN", "amine", 1, "ethanamine", "乙胺"),
    ("CC(N)C", "amine", 1, "propan-2-amine", "丙-2-胺"),
    ("NCCN", "amine", 2, "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("NCCCN", "amine", 2, "propane-1,3-diamine", "丙烷-1,3-二胺"),
    ("NCC(N)CN", "amine", 3, "propane-1,2,3-triamine", "丙烷-1,2,3-三胺"),
]


def p62_primary_polyamine_principal_expression_facts___typed_parent(smiles, kind, count):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("kind") == kind
                and p.get("principal_expression_facts")
                and p["principal_expression_facts"].multiplicity == count)


@pytest.mark.parametrize("smiles,kind,count,en,zh", p62_primary_polyamine_principal_expression_facts__CASES)
def test_primary_polyamine_typed_facts_and_names(smiles, kind, count, en, zh):
    info = analyze(Chem.MolFromSmiles(smiles))
    facts = p62_primary_polyamine_principal_expression_facts___typed_parent(smiles, kind, count)["principal_expression_facts"]
    assert facts.group_class.value == "amine"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert len(facts.attachment_atoms) == count
    expected = {a for o in inventory_from_info(info).occurrences(FG.AMINE) for a in o.parent_anchors}
    assert set(facts.attachment_atoms) == expected
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


# ==========================================================================
# 合并自 test_p63_alcohol_principal_expression_facts.py
# IUPAC: P-63.1.1 / P-63.1.2
# Layer: L2,L4,L5
# ==========================================================================
p63_alcohol_principal_expression_facts__CASES = [
    ("CCCO", "alcohol", 1, "propan-1-ol", "丙-1-醇"),
    ("CC(O)C", "alcohol", 1, "propan-2-ol", "丙-2-醇"),
    ("OCCCO", "alcohol", 2, "propane-1,3-diol", "丙烷-1,3-二醇"),
    ("CC(O)C(O)C", "alcohol", 2, "butane-2,3-diol", "丁烷-2,3-二醇"),
    ("OCC(O)CO", "alcohol", 3, "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("OC1CCC(O)CC1", "alcohol", 2, "cyclohexane-1,4-diol", "环己烷-1,4-二醇"),
]


def p63_alcohol_principal_expression_facts___typed_parent(smiles: str, kind: str, count: int):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("kind") == kind
                and p.get("principal_expression_facts")
                and p["principal_expression_facts"].multiplicity == count)


@pytest.mark.parametrize("smiles,kind,count,en,zh", p63_alcohol_principal_expression_facts__CASES)
def test_alcohol_principal_expression_facts_and_names(smiles, kind, count, en, zh):
    parent = p63_alcohol_principal_expression_facts___typed_parent(smiles, kind, count)
    facts = parent["principal_expression_facts"]
    assert facts.group_class.value == "alcohol"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert len(facts.attachment_atoms) == count
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


# ==========================================================================
# 合并自 test_p63_ketone_principal_expression_facts.py
# IUPAC: P-63.2.1
# Layer: L2,L4,L5
# ==========================================================================
p63_ketone_principal_expression_facts__CASES = [
    ("CC(=O)CC", "ketone", 1, "butan-2-one", "丁-2-酮"),
    ("CCC(=O)CC", "ketone", 1, "pentan-3-one", "戊-3-酮"),
    ("CC(=O)CC(=O)C", "ketone", 2, "pentane-2,4-dione", "戊-2,4-二酮"),
    ("CC(=O)C(C)=O", "ketone", 2, "butane-2,3-dione", "丁-2,3-二酮"),
]


def p63_ketone_principal_expression_facts___parent(smiles, kind):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(parent for parent in parents if parent.get("kind") == kind)


@pytest.mark.parametrize("smiles,kind,count,en,zh", p63_ketone_principal_expression_facts__CASES)
def test_ketone_principal_expression_facts_and_names(smiles, kind, count, en, zh):
    facts = p63_ketone_principal_expression_facts___parent(smiles, kind)["principal_expression_facts"]
    assert facts.group_class.value == "ketone"
    assert facts.multiplicity == count
    assert facts.relation is PrincipalRelation.IN_SKELETON
    assert len(facts.occurrence_ids) == count
    assert facts.characteristic_atoms and facts.attachment_atoms
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_alcohol_near_neighbor_has_no_ketone_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCCCO")))
    alcohol = next(parent for parent in parents if parent.get("kind") == "alcohol")
    facts = alcohol.get("principal_expression_facts")
    assert facts is None or facts.group_class.value != "ketone"
