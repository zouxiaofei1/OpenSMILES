# IUPAC: P-44
# Layer: L1,L2
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer


CASES = [
    # 苯环 + FG 的 kind 收敛为 FG 类别；保留名由 L5 typed_kinds 决定。
    ("O=C(O)c1ccccc1", "acid", "benzoic acid"),
    ("O=Cc1ccccc1", "aldehyde", "benzaldehyde"),
    ("N#Cc1ccccc1", "nitrile", "benzonitrile"),
    ("NC(=O)c1ccccc1", "amide", "benzamide"),
    ("Oc1ccccc1", "alcohol", "phenol"),
    ("Nc1ccccc1", "amine", "aniline"),
    ("CCCC(=O)OC", "ester", "methyl butanoate"),
    ("CCCC(=O)N", "amide", "butanamide"),
    ("CCCC=O", "aldehyde", "butanal"),
    ("CC(=O)CC", "ketone", "butan-2-one"),
    ("CCCC#N", "nitrile", "butanenitrile"),
]


@pytest.mark.parametrize("smiles,kind,en", CASES)
def test_topology_first_principal_path_owns_target(smiles, kind, en):
    # fg_try_fns 已随 fg 注册层删除；principal-only 是默认状态。
    info = analyze(Chem.MolFromSmiles(smiles))
    parents = _collect_candidates(info)
    result = SMILESNNamer().name(smiles)
    owned = [parent for parent in parents if parent.get("kind") == kind]
    assert owned
    count = 2 if kind == "diester" else 1
    assert owned[0].get("principal_group_count") == count
    assert len(owned[0].get("covered_principal_ids") or ()) == count
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


def test_open_chain_monoester_is_not_claimed_as_diester():
    info = analyze(Chem.MolFromSmiles("COC(=O)CCC"))
    assert all(parent.get("kind") != "diester" for parent in _collect_candidates(info))


@pytest.mark.parametrize("smiles,kind,group_class", [
    # 环 + FG 一律收敛为 FG 类别 kind（苯/饱和环/稠环/杂环平等，词干由 scaffold 承载）。
    ("O=C(O)C1CCCCC1", "acid", "acid"),
    ("N#CC1CCCCC1", "nitrile", "nitrile"),
    ("O=C(O)c1ccc2ccccc2c1", "acid", "acid"),
])
def test_ring_principal_uses_base_scaffold_kind(smiles, kind, group_class):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    matched = [parent for parent in parents if parent.get("kind") == kind]
    assert matched
    assert matched[0]["principal_group_class"] == group_class
    assert "carboxylic" not in matched[0]["kind"]
    assert "carbonitrile" not in matched[0]["kind"]
