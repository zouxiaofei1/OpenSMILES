"""Unit tests for alkenyl topology probes (L2 side_alkyl)."""
from __future__ import annotations

from rdkit.Chem import MolFromSmiles, GetSymmSSSR

from namepredict.tools.side_alkyl import (
    _is_vinyl,
    _is_allyl,
    _is_isopropenyl,
)


def _ring_set(mol):
    rings = GetSymmSSSR(mol)
    return set(rings[0]) if rings else set()


# ── vinyl (Parent–CH=CH2) ──

def test_is_vinyl_detects_ethenyl_side_chain():
    """C=CC1CCCCC1: cyclohexane-vinyl."""
    mol = MolFromSmiles("C=CC1CCCCC1")
    ring = _ring_set(mol)
    # vinyl start is atom 1 (=CH– ring)
    atoms = _is_vinyl(mol, 1, ring)
    assert atoms is not None
    assert len(atoms) == 2


def test_is_vinyl_returns_none_for_saturated():
    mol = MolFromSmiles("CCC1CCCCC1")
    atoms = _is_vinyl(mol, 1, _ring_set(mol))
    assert atoms is None


# ── allyl (Parent–CH2–CH=CH2) ──

def test_is_allyl_detects_propenyl_side_chain():
    """C=CCC1CCCCC1: cyclohexane-allyl."""
    mol = MolFromSmiles("C=CCC1CCCCC1")
    ring = _ring_set(mol)
    atoms = _is_allyl(mol, 2, ring)
    assert atoms is not None
    assert len(atoms) == 3


def test_is_allyl_returns_none_for_saturated():
    mol = MolFromSmiles("CCCC1CCCCC1")
    atoms = _is_allyl(mol, 1, _ring_set(mol))
    assert atoms is None


# ── isopropenyl (Parent–C(CH3)=CH2) ──

def test_is_isopropenyl_detects_methylvinyl_side_chain():
    """CC(=C)C1CCCCC1: cyclohexane-isopropenyl."""
    mol = MolFromSmiles("CC(=C)C1CCCCC1")
    ring = _ring_set(mol)
    atoms = _is_isopropenyl(mol, 1, ring)
    assert atoms is not None
    assert len(atoms) == 3


def test_is_isopropenyl_returns_none_for_n_propyl():
    mol = MolFromSmiles("CCCC1CCCCC1")
    atoms = _is_isopropenyl(mol, 1, _ring_set(mol))
    assert atoms is None
