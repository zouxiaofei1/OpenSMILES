# IUPAC: P-62.2.2.1
# Layer: L2,L3,L4,L5
import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import _collect_candidates
from namepredict.namer import SMILESNNamer


def _arm_parent(smiles):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    return next(p for p in parents if p.get("neutral_amine_arm_facts"))


def test_primary_amine_has_no_n_arm_facts():
    parents = _collect_candidates(analyze(Chem.MolFromSmiles("CCN")))
    assert not any(p.get("neutral_amine_arm_facts") for p in parents)
