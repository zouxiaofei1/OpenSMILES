# IUPAC: P-24.2.1
# Layer: L2
"""Spiro parent hydride detection — Round 1: saturated all-carbon monospiro."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.spiro_parent import _is_simple_spiro, _try_spiro_parent


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


# ── Positive cases ──


def test_spiro45_decane_detected():
    info = _info("C1CCC2(C1)CCCCC2")
    sp = _is_simple_spiro(info)
    assert sp is not None
    assert sp["topology"] == "spiro"
    assert sp["ring_sizes"] == [4, 5]


def test_spiro44_nonane_detected():
    info = _info("C1CCC2(CCCC2)C1")
    sp = _is_simple_spiro(info)
    assert sp is not None
    assert sp["topology"] == "spiro"
    assert sp["ring_sizes"] == [4, 4]


def test_spiro55_undecane_detected():
    info = _info("C1CCCC2(CCCCC2)C1")
    sp = _is_simple_spiro(info)
    assert sp is not None
    assert sp["topology"] == "spiro"
    assert sp["ring_sizes"] == [5, 5]


def test_try_spiro_returns_parent():
    info = _info("C1CCC2(C1)CCCCC2")
    parent = _try_spiro_parent(info)
    assert parent is not None
    assert parent["kind"] == "spiro"
    assert parent["stem_en"] == "spiro[4.5]"
    assert parent["stem_zh"] == "螺[4.5]"
    assert parent["ring_sizes"] == [4, 5]


# ── Negative cases ──


def test_benzene_not_spiro():
    info = _info("c1ccccc1")
    assert _is_simple_spiro(info) is None


def test_fused_not_spiro():
    info = _info("c1ccc2ccccc2c1")  # naphthalene
    assert _is_simple_spiro(info) is None


def test_open_chain_not_spiro():
    info = _info("CCCCCC")
    assert _is_simple_spiro(info) is None


def test_biphenyl_not_spiro():
    info = _info("c1ccc(-c2ccccc2)cc1")
    assert _is_simple_spiro(info) is None


@pytest.mark.parametrize("smiles", [
    "c1ccc2ccccc2c1",       # fused naphthalene
    "CCCCCC",                # open chain
    "c1ccccc1",              # benzene
    "C1CCCCC1",              # cyclohexane
    "C1CCOC1",               # oxolane (hetero)
])
def test_non_spiro_returns_none(smiles):
    assert _try_spiro_parent(_info(smiles)) is None
