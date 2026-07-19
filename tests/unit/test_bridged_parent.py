# IUPAC: P-23 (von Baeyer bridged parent hydrides)
# Layer: L1, L2
"""Bridged (von Baeyer) parent hydride detection — Round 1: saturated all-carbon bicyclic."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.bridged_parent import _is_simple_bridged, _try_bridged_parent


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


# ── Positive cases ──


def test_norbornane_detected():
    """bicyclo[2.2.1]heptane"""
    info = _info("C1CC2CCC1C2")
    br = _is_simple_bridged(info)
    assert br is not None
    assert br["topology"] == "bridged"
    assert br["bridgeheads"] == [2, 5]
    assert br["bridge_lengths"] == [2, 2, 1]


def test_bicyclo321_octane_detected():
    """bicyclo[3.2.1]octane"""
    info = _info("C1CC2CCC(C1)C2")
    br = _is_simple_bridged(info)
    assert br is not None
    assert br["topology"] == "bridged"
    assert br["bridge_lengths"] == [3, 2, 1]


def test_try_bridged_returns_parent():
    info = _info("C1CC2CCC1C2")
    parent = _try_bridged_parent(info)
    assert parent is not None
    assert parent["kind"] == "bridged"
    assert parent["stem_en"] == "bicyclo[2.2.1]"
    assert parent["stem_zh"] == "双环[2.2.1]"
    assert parent["bridge_lengths"] == [2, 2, 1]
    assert parent["bridgeheads"] == [2, 5]


def test_bicyclo321_stem():
    info = _info("C1CC2CCC(C1)C2")
    parent = _try_bridged_parent(info)
    assert parent is not None
    assert parent["stem_en"] == "bicyclo[3.2.1]"
    assert parent["bridge_lengths"] == [3, 2, 1]


# ── Negative cases ──


def test_benzene_not_bridged():
    info = _info("c1ccccc1")
    assert _is_simple_bridged(info) is None


def test_fused_not_bridged():
    info = _info("c1ccc2ccccc2c1")  # naphthalene
    assert _is_simple_bridged(info) is None


def test_spiro_not_bridged():
    info = _info("C1CCC2(C1)CCCCC2")  # spiro[4.5]decane
    assert _is_simple_bridged(info) is None


def test_open_chain_not_bridged():
    info = _info("CCCCCC")
    assert _is_simple_bridged(info) is None


def test_cyclohexane_not_bridged():
    info = _info("C1CCCCC1")
    assert _is_simple_bridged(info) is None


@pytest.mark.parametrize("smiles", [
    "c1ccc2ccccc2c1",       # fused naphthalene
    "CCCCCC",                # open chain
    "c1ccccc1",              # benzene
    "C1CCCCC1",              # cyclohexane
    "C1CCOC1",               # oxolane (hetero)
    "C1CCC2(C1)CCCCC2",     # spiro[4.5]decane
])
def test_non_bridged_returns_none(smiles):
    assert _try_bridged_parent(_info(smiles)) is None


# ── Integration: full pipeline ──


def test_norbornane_full_pipeline():
    """End-to-end: norbornane → bicyclo[2.2.1]heptane / 双环[2.2.1]庚烷"""
    from namepredict.layer2.parent_selector import select_parent
    from namepredict.layer4.numbering import number
    from namepredict.layer5.assembler import assemble

    info = _info("C1CC2CCC1C2")
    parent = select_parent(info)
    assert parent["kind"] == "bridged"
    numbered = number(parent, info.get("substituents", []))
    result = assemble(numbered)
    assert result.success
    assert result.en == "bicyclo[2.2.1]heptane"
    assert result.zh == "双环[2.2.1]庚烷"


def test_norbornane_numbering_chain():
    """von Baeyer numbering: 7 unique atoms in correct order."""
    from namepredict.layer2.parent_selector import select_parent
    from namepredict.layer4.numbering import number

    info = _info("C1CC2CCC1C2")
    parent = select_parent(info)
    numbered = number(parent, info.get("substituents", []))
    chain = numbered["parent"]["chain"]
    assert len(chain) == 7
    assert len(set(chain)) == 7  # all unique
    # First two positions should be bridgeheads
    assert chain[0] in parent["bridgeheads"]
    assert chain[3] in parent["bridgeheads"]  # second bridgehead at position 4


def test_bicyclo321_full_pipeline():
    """End-to-end: bicyclo[3.2.1]octane"""
    from namepredict.layer2.parent_selector import select_parent
    from namepredict.layer4.numbering import number
    from namepredict.layer5.assembler import assemble

    info = _info("C1CC2CCC(C1)C2")
    parent = select_parent(info)
    assert parent["kind"] == "bridged"
    numbered = number(parent, info.get("substituents", []))
    result = assemble(numbered)
    assert result.success
    assert result.en == "bicyclo[3.2.1]octane"
    assert "双环[3.2.1]" in result.zh
