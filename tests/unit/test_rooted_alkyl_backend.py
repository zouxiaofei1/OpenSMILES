"""Rooted-tree backend: pure saturated carbon claims only."""
from __future__ import annotations

from rdkit import Chem

from namepredict.layer1.analyzer import analyze
from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot, iter_claims
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer3.substituent_namer import RootedTreeBackend, SubstituentNamer


def _claim(slot: SideSlot, attach: int, root: int, atoms: set[int]) -> ClaimedBlock:
    return ClaimedBlock(
        slot=slot,
        attach_parent=attach,
        root=root,
        atoms=frozenset(atoms),
    )


def _benzene_methylbutyl_claim():
    """(2-methylbutyl)benzene side: ring owned, five side carbons claimed."""
    mol = Chem.MolFromSmiles("c1ccc(cc1)CC(C)CC")
    # ring atoms 0-5; side: 6(CH2 root)-7(CH)-8(Me)-9(CH2)-10(Me)
    claim = _claim(SideSlot.RING_C, attach=3, root=6, atoms={6, 7, 8, 9, 10})
    return mol, claim


def _acid_methylbutyl_claim():
    mol = Chem.MolFromSmiles("CCCCC(CC(C)CC)C(=O)O")
    info = analyze(mol)
    parent = select_parent(info)
    claims = iter_claims(mol, parent["owned_atoms"])
    assert len(claims) == 1
    return mol, claims[0]


def test_rooted_tree_names_2_methylbutyl_on_benzene():
    mol, claim = _benzene_methylbutyl_claim()
    hit = RootedTreeBackend().try_name(mol, claim, depth=0)
    assert hit is not None
    assert hit.backend == "rooted_tree"
    assert hit.claim is claim
    assert hit.en == "2-methylbutyl"
    assert hit.zh == "2-甲基丁基"


def test_rooted_tree_names_2_methylbutyl_on_acid():
    mol, claim = _acid_methylbutyl_claim()
    hit = RootedTreeBackend().try_name(mol, claim, depth=0)
    assert hit is not None
    assert hit.backend == "rooted_tree"
    assert "2-methylbutyl" in hit.en
    assert "2-甲基丁基" in hit.zh


def test_alkene_side_returns_none_from_rooted_tree():
    """Allyl (unsaturated) must not be accepted by rooted-tree backend."""
    mol = Chem.MolFromSmiles("C=CCc1ccccc1")
    # force a claim over the allyl carbons only (0,1,2)
    claim = _claim(SideSlot.RING_C, attach=3, root=2, atoms={0, 1, 2})
    assert RootedTreeBackend().try_name(mol, claim, depth=0) is None


def test_over_limit_13_atom_tree_returns_none():
    """13 saturated carbons from root exceeds max_atoms=12."""
    # Ph-(CH2)12-CH3 → 13 side carbons
    mol = Chem.MolFromSmiles("c1ccccc1CCCCCCCCCCCCC")
    side = {i for i in range(6, 19)}  # 13 carbons
    claim = _claim(SideSlot.RING_C, attach=0, root=6, atoms=side)
    assert len(claim.atoms) == 13
    assert RootedTreeBackend().try_name(mol, claim, depth=0) is None


def test_namer_order_reaches_rooted_tree_before_recursive():
    mol, claim = _benzene_methylbutyl_claim()
    hit = SubstituentNamer().name(mol, claim, depth=0)
    assert hit is not None
    assert hit.backend == "rooted_tree"
    assert hit.en == "2-methylbutyl"
