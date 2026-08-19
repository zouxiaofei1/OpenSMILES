# IUPAC: P-25.3.3
# Layer: L4
"""fused_numbering: 外周骨架编号 + 稠合碳 a/b/c 字母位次(P-25.3.3.1)。"""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer4.fused_orientation import preferred_orientation
from namepredict.layer4.fused_numbering import number_fused_system
from namepredict.layer2.ring_scaffold import FUSED56_LABELS, NAPH_LABELS


def _number(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    systems = build_ring_systems(mol)
    assert systems
    rings = list(mol.GetRingInfo().AtomRings())
    orient = preferred_orientation(mol, rings, systems[0]["fusion_edges"])
    assert orient is not None
    return mol, number_fused_system(mol, rings, orient.coord_dict())


def test_naphthalene_labels_match_standard():
    _, (chain, labels) = _number("c1ccc2ccccc2c1")
    assert len(chain) == 10
    assert labels == list(NAPH_LABELS)


def test_indole_labels_match_standard():
    _, (chain, labels) = _number("c1ccc2[nH]ccc2c1")
    assert len(chain) == 9
    assert labels == list(FUSED56_LABELS)


def test_quinoline_hetero_atom_locant_one():
    mol, (chain, labels) = _number("c1ccc2ncccc2c1")
    n = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() == 7][0]
    assert labels[chain.index(n)] == "1"


def test_phenanthrene_14_atoms_with_letter_fused_carbons():
    mol, (chain, labels) = _number("c1ccc2c(c1)ccc1ccccc12")
    assert len(chain) == 14
    fused = {a for a in chain if labels[chain.index(a)][-1].isalpha()}
    assert len(fused) >= 3  # 稠合碳均有字母位次


def test_azulene_odd_ring_pair_numbering():
    """azulene 5+7 奇环对: 10 原子全部编号, 桥头字母位。"""
    _, (chain, labels) = _number("C1=CC=C2C=CC=CC2=C1")
    assert len(chain) == 10
    assert any(lbl[-1].isalpha() for lbl in labels)
