# IUPAC: P-41 cation / tetraalkylammonium
# Layer: L1,L2,L5
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer


def test_quaternary_is_not_neutral_amine():
    info = analyze(Chem.MolFromSmiles("C[N+](C)(C)C"))
    assert not info["amines"]
