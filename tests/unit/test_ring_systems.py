# IUPAC: P-25 / P-22
# Layer: L1
"""Ring-system topology: fusion components for naphthalene/indole/quinoline/spiro."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze


def _systems(smiles: str) -> list[dict]:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)["ring_systems"]


def test_benzene_mono():
    syss = _systems("c1ccccc1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 1
    assert s["topology"] == "mono"
    assert s["n_atoms"] == 6
    assert s["is_aromatic_mancude"] is True
    assert s["hetero_atoms"] == []


def test_naphthalene_fused2():
    syss = _systems("c1ccc2ccccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["topology"] == "fused"
    assert s["n_atoms"] == 10
    assert len(s["fusion_edges"]) == 1
    assert len(s["fusion_edges"][0][2]) == 2  # shared atoms
    assert s["is_aromatic_mancude"] is True


def test_indole_fused56():
    syss = _systems("c1ccc2[nH]ccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["topology"] == "fused"
    assert s["n_atoms"] == 9
    zs = {h["Z"] for h in s["hetero_atoms"]}
    assert 7 in zs


def test_quinoline_fused66():
    syss = _systems("c1ccc2ncccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["topology"] == "fused"
    assert s["n_atoms"] == 10
    assert any(h["Z"] == 7 for h in s["hetero_atoms"])


def test_biphenyl_two_systems():
    """Two unfused rings → two mono systems (not fused)."""
    syss = _systems("c1ccc(-c2ccccc2)cc1")
    assert len(syss) == 2
    assert all(s["topology"] == "mono" for s in syss)
    assert all(s["n_rings"] == 1 for s in syss)


def test_spiro45_marked():
    """Spiro[4.5]decane: two rings share 1 atom → separate mono systems."""
    syss = _systems("C1CCC2(C1)CCCCC2")
    assert len(syss) == 2
    assert all(s["n_rings"] == 1 for s in syss)
    # spiro pairs exist but each system is mono (not fused)
    assert all(s["topology"] in ("mono", "spiro", "other") for s in syss)


def test_pyridine_hetero():
    syss = _systems("c1ccncc1")
    assert len(syss) == 1
    s = syss[0]
    assert s["topology"] == "mono"
    assert any(h["Z"] == 7 for h in s["hetero_atoms"])


def test_open_chain_empty():
    assert _systems("CCCC") == []
