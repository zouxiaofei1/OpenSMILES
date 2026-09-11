# IUPAC: P-25.3.3
# Layer: L4
"""fused_numbering: 外周骨架编号 + 稠合碳 a/b/c 字母位次(P-25.3.3.1)。"""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer4.fused_orientation import preferred_orientations
from namepredict.layer4.fused_numbering import fused_atoms, number_fused_system
from namepredict.layer2.ring_scaffold import FUSED56_LABELS, NAPH_LABELS


def _number(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    systems = build_ring_systems(mol)
    assert systems
    rings = list(mol.GetRingInfo().AtomRings())
    orients = preferred_orientations(mol, rings, systems[0]["fusion_edges"])
    assert orients
    return mol, number_fused_system(mol, rings, [o.coord_dict() for o in orients])


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


def test_walk_starts_next_to_a_fused_atom():
    """P-25.3.3.1.1 外周行走自稠合边一端起算: 起点取偏环顶端的顶点会使位次整体错一位。"""
    for smiles in ("COc1cc2oc(=O)c3c(O)cc(O)cc3c2c(C)c1Cl",
                   "C1=C2C=3N(C(NC2=CC=C1)=O)C1=C(N3)C=CC=C1"):
        mol, (chain, _) = _number(smiles)
        fused = fused_atoms(list(mol.GetRingInfo().AtomRings()))
        ring = next(r for r in mol.GetRingInfo().AtomRings() if chain[0] in r)
        i = ring.index(chain[0])
        assert ring[i - 1] in fused or ring[(i + 1) % len(ring)] in fused


def test_benzo_c_chromenone_lactone_locant():
    """苯并[c]色烯-6-酮的环内羰基得 6 位(起点紧邻稠合边; 起点偏顶点时错编为 5-酮)。"""
    mol, (chain, labels) = _number("COc1cc2oc(=O)c3c(O)cc(O)cc3c2c(C)c1Cl")
    locants = dict(zip(chain, labels))
    carbonyl = next(a for a in chain if any(
        n.GetSymbol() == "O" and mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondType().name == "DOUBLE"
        for n in mol.GetAtomWithIdx(a).GetNeighbors()))
    assert locants[carbonyl] == "6"


def test_6_5_6_three_ring_numbering():
    """6-5-6 线性稠环(变形五元环居中): 三个 N 得最低位次, 稠合碳得字母位次(a/b/c)。"""
    mol, (chain, labels) = _number("ClC1C2NC3C=CC=CC=3C=2N=CN=1")
    assert len(chain) == 13
    heteros = sorted(labels[chain.index(a)] for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
    assert heteros == ["1", "3", "5"]  # P-25.3.3.2.3 杂原子最低位次
    fused_letter = [lbl for lbl in labels if lbl[-1].isalpha()]
    assert len(fused_letter) >= 4  # 稠合碳均有字母位次
