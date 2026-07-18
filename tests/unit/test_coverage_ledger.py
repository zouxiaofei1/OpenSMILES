"""Coverage ledger: heavy-atom gap/overlap over ownership + named claims."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer3.coverage import CoverageLedger, build_coverage_ledger
from namepredict.layer3.substituent_namer import SubstituentName


def _heavy(mol) -> frozenset[int]:
    return frozenset(
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1
    )


def _name(claim: ClaimedBlock, en: str = "x") -> SubstituentName:
    return SubstituentName(
        claim=claim,
        en=en,
        zh=en,
        requires_parentheses=False,
        backend="retained",
    )


def _claim(atoms: frozenset[int], attach: int = 0, root: int = 1) -> ClaimedBlock:
    return ClaimedBlock(
        slot=SideSlot.CHAIN_C,
        attach_parent=attach,
        root=root,
        atoms=atoms,
    )


def test_complete_ledger_when_owned_and_names_cover_all_heavy():
    """Disjoint owned + named claims covering all heavy atoms → complete."""
    mol = preprocess("CC")  # two carbons
    owned = frozenset([0])
    names = [_name(_claim(frozenset([1]), attach=0, root=1))]
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=names)
    assert isinstance(ledger, CoverageLedger)
    assert ledger.owned_atoms == owned
    assert ledger.named_claims == tuple(names)
    assert ledger.gap == frozenset()
    assert ledger.overlap == frozenset()
    assert ledger.complete is True
    assert _heavy(mol) == frozenset([0, 1])


def test_gap_when_one_heavy_atom_omitted():
    """Omitting one heavy atom places it in gap; not complete."""
    mol = preprocess("CCC")  # three carbons
    owned = frozenset([0])
    names = [_name(_claim(frozenset([1]), attach=0, root=1))]
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=names)
    assert ledger.gap == frozenset([2])
    assert ledger.overlap == frozenset()
    assert ledger.complete is False


def test_overlap_when_atom_in_owned_and_claim():
    """Atom in both owned_atoms and a named claim → overlap."""
    mol = preprocess("CC")
    owned = frozenset([0, 1])
    names = [_name(_claim(frozenset([1]), attach=0, root=1))]
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=names)
    assert ledger.overlap == frozenset([1])
    assert ledger.gap == frozenset()
    assert ledger.complete is False


def test_overlap_when_atom_in_two_named_claims():
    """Atom appearing in two named claims → overlap."""
    mol = preprocess("CCC")
    owned = frozenset([0])
    names = [
        _name(_claim(frozenset([1, 2]), attach=0, root=1), en="a"),
        _name(_claim(frozenset([2]), attach=0, root=2), en="b"),
    ]
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=names)
    assert 2 in ledger.overlap
    assert ledger.complete is False


def test_hydrogen_excluded_from_gap():
    """Hydrogens are never counted in gap or overlap."""
    mol = preprocess("C")  # methane: 1 C + 4 H
    owned = frozenset([0])
    ledger = build_coverage_ledger(mol, owned_atoms=owned, names=[])
    assert ledger.gap == frozenset()
    assert ledger.overlap == frozenset()
    assert ledger.complete is True
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() != 1 for i in ledger.gap)
