# IUPAC: P-62.2.2.1
# Layer: L2,L3,L4,L5
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer


CASES = [
    ("CCNCC", 2, 1, "N-ethylethanamine", "N-乙基乙胺"),
    ("CCNC", 2, 1, "N-methylethanamine", "N-甲基乙胺"),
    ("CCN(CC)CC", 3, 2, "N,N-diethylethanamine", "N,N-二乙基乙胺"),
    ("CN(C)C", 3, 2, "N,N-dimethylmethanamine", "N,N-二甲基甲胺"),
]


def _arm_parent(smiles):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("neutral_amine_arm_facts"))


@pytest.mark.parametrize("smiles,degree,n_subs,en,zh", CASES)
def test_neutral_amine_arm_facts_and_names(smiles, degree, n_subs, en, zh):
    parent = _arm_parent(smiles)
    facts = parent["neutral_amine_arm_facts"]
    assert facts.degree == degree
    assert len(facts.substituent_arms) == n_subs
    assert tuple(parent["chain"]) == facts.parent_arm.atom_ids
    assert all(arm.atom_ids for arm in facts.substituent_arms)
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_primary_amine_has_no_n_arm_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCN")))
    assert not any(p.get("neutral_amine_arm_facts") for p in parents)
