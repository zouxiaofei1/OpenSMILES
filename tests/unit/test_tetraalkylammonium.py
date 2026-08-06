# IUPAC: P-41 cation / tetraalkylammonium
# Layer: L1,L2,L5
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer

CASES = [
    ("C[N+](C)(C)C", "tetramethylammonium", "四甲基铵"),
    ("CC[N+](C)(C)C", "ethyltrimethylammonium", "乙基三甲基铵"),
    ("CC[N+](CC)(C)C", "diethyldimethylammonium", "二乙基二甲基铵"),
]

@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_tetraalkylammonium_names(smiles, en, zh):
    info = analyze(Chem.MolFromSmiles(smiles))
    assert len(info["quaternary_ammoniums"]) == 1
    parent = next(p for p in _collect_candidates(info) if p.get("kind") == "tetraalkylammonium")
    assert len(parent["quaternary_arm_lengths"]) == 4
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == en
    assert normalize_zh(result.zh) == zh


def test_quaternary_is_not_neutral_amine():
    info = analyze(Chem.MolFromSmiles("C[N+](C)(C)C"))
    assert not info["amines"]
