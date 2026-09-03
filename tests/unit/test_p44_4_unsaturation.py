# IUPAC: P-44.4.1.1 / P-44.4.1.2
# Layer: L2
from rdkit import Chem

from namepredict.layer2.parent_skeleton import (
    ParentSkeleton,
    SkeletonTopology,
    keep_p44_4_unsaturation,
)


def _chain(atoms):
    return ParentSkeleton(SkeletonTopology.ACYCLIC, tuple(atoms), frozenset({"fg:0"}))


def test_p44_4_prefers_more_multiple_bonds():
    mol = Chem.MolFromSmiles("C=CC.C=CC#C")
    one, two = _chain(range(3)), _chain(range(3, 7))
    assert keep_p44_4_unsaturation(mol, (one, two)) == (two,)


def test_p44_4_prefers_more_double_bonds_when_total_ties():
    mol = Chem.MolFromSmiles("C=CC=C.C=CC#C")
    diene, enyne = _chain(range(4)), _chain(range(4, 8))
    assert keep_p44_4_unsaturation(mol, (enyne, diene)) == (diene,)


def test_p44_4_aromatic_counts_as_kekule_multiple_bonds():
    mol = Chem.MolFromSmiles("c1ccccc1.C1CCCCC1")
    aromatic, saturated = _chain(range(6)), _chain(range(6, 12))
    assert keep_p44_4_unsaturation(mol, (aromatic, saturated)) == (aromatic,)


def test_p44_4_pyrimidine_beats_saturated_ring():
    # 哌嗪-嘧啶环环相连：嘧啶(芳香,3 多重键当量)应优先于哌嗪(饱和,0)。
    mol = Chem.MolFromSmiles("C1C(N2C(C)CNCC2)=NC=NC=1")
    pymidine = ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple([0, 1, 9, 10, 11, 12]), frozenset())
    piperazine = ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple([2, 3, 5, 6, 7, 8]), frozenset())
    assert keep_p44_4_unsaturation(mol, (pymidine, piperazine)) == (pymidine,)
