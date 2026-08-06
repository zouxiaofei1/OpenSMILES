# IUPAC: P-44.4.1.1 / P-44.4.1.2
# Layer: L2
from rdkit import Chem

from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, keep_p44_4_unsaturation


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


def test_p44_4_ignores_aromatic_bonds_and_preserves_ties():
    mol = Chem.MolFromSmiles("c1ccccc1.C1CCCCC1")
    aromatic, saturated = _chain(range(6)), _chain(range(6, 12))
    assert keep_p44_4_unsaturation(mol, (aromatic, saturated)) == (aromatic, saturated)
