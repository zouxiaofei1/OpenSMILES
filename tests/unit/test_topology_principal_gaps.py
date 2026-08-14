# IUPAC: P-44
# Layer: L1,L2
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer


CASES = [
    ("O=C(O)c1ccccc1", "benzoic", "benzoic acid"),
    ("O=Cc1ccccc1", "benzaldehyde", "benzaldehyde"),
    ("N#Cc1ccccc1", "benzonitrile", "benzonitrile"),
    ("NC(=O)c1ccccc1", "benzamide", "benzamide"),
    # phenol 保留名迁 L5：kind=benzene，输出仍 phenol
    ("Oc1ccccc1", "benzene", "phenol"),
    ("Nc1ccccc1", "aniline", "aniline"),
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
    ("O=C(O)C1CCCCC1", "cycloalkane", "acid"),
    ("N#CC1CCCCC1", "cycloalkane", "nitrile"),
    ("O=C(O)c1ccc2ccccc2c1", "naphthalene", "acid"),
])
def test_ring_principal_uses_base_scaffold_kind(smiles, kind, group_class):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    matched = [parent for parent in parents if parent.get("kind") == kind]
    assert matched
    assert matched[0]["principal_group_class"] == group_class
    assert "carboxylic" not in matched[0]["kind"]
    assert "carbonitrile" not in matched[0]["kind"]
