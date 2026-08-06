# IUPAC: P-44.2
# Layer: L2
from rdkit import Chem

from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, keep_p44_2


def _ring(atoms):
    return ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(atoms), frozenset({"fg:0"}))


def test_p44_2_prefers_heterocycle_over_carbocycle():
    mol = Chem.MolFromSmiles("C1CC1.N1CC1")
    assert keep_p44_2(mol, (_ring((0, 1, 2)), _ring((3, 4, 5)))) == (_ring((3, 4, 5)),)


def test_p44_2_prefers_more_nitrogens():
    mol = Chem.MolFromSmiles("N1CCCC1.N1CCNC1")
    one = _ring(range(5))
    two = _ring(range(5, 10))
    assert keep_p44_2(mol, (one, two)) == (two,)


def test_p44_2_without_n_prefers_senior_heteroatom():
    mol = Chem.MolFromSmiles("O1CCC1.S1CCC1")
    oxygen = _ring(range(4))
    sulfur = _ring(range(4, 8))
    assert keep_p44_2(mol, (sulfur, oxygen)) == (oxygen,)


def test_p44_2_prefers_more_rings_before_more_atoms():
    mol = Chem.MolFromSmiles("N1CCC2CC12.N1CCCCC1")
    bicyclic = _ring(range(6))
    monocyclic = _ring(range(6, 12))
    assert keep_p44_2(mol, (monocyclic, bicyclic)) == (bicyclic,)


def test_p44_2_prefers_more_ring_atoms_after_hetero_ties():
    mol = Chem.MolFromSmiles("N1CCC1.N1CCCC1")
    small = _ring(range(4))
    large = _ring(range(4, 9))
    assert keep_p44_2(mol, (small, large)) == (large,)
