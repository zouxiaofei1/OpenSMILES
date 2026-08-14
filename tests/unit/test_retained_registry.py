# IUPAC: P-22 / P-25
# Layer: L2
"""Retained registry matches ring_systems for benzene/pyridine/naphthalene/indole."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.ring_scaffold import match_scaffold_ids, match_systems


def _ids(smiles: str) -> list[str]:
    mol = preprocess(smiles)
    assert mol is not None
    return match_scaffold_ids(analyze(mol))


def test_benzene_registry():
    assert "benzene" in _ids("c1ccccc1")


def test_pyridine_registry():
    assert "pyridine" in _ids("c1ccncc1")


def test_naphthalene_registry():
    assert "naphthalene" in _ids("c1ccc2ccccc2c1")


def test_indole_registry():
    assert "indole" in _ids("c1ccc2[nH]ccc2c1")


def test_biphenyl_two_benzene():
    ids = _ids("c1ccc(-c2ccccc2)cc1")
    assert ids.count("benzene") == 2


def test_open_chain_no_match():
    assert _ids("CCCC") == []


def test_match_systems_has_entry():
    mol = preprocess("c1ccc2ccccc2c1")
    info = analyze(mol)
    hits = match_systems(info)
    assert len(hits) == 1
    sid, system, entry = hits[0]
    assert sid == "naphthalene"
    assert entry["en"] == "naphthalene"
    assert system["n_rings"] == 2
