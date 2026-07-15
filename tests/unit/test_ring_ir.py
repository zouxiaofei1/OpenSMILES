# IUPAC: P-25 / P-22
# Layer: L1
"""RingSystemIR: typed SSSR fusion topology (benzene / naphthalene / indole)."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_ir import RingSystemIR, build_ring_ir


def _ir(smiles: str) -> list[RingSystemIR]:
    mol = preprocess(smiles)
    assert mol is not None
    return build_ring_ir(mol)


def test_benzene_mono():
    syss = _ir("c1ccccc1")
    assert len(syss) == 1
    s = syss[0]
    assert len(s.components) == 1
    assert s.topology == "mono"
    assert len(s.atom_ids) == 6
    assert s.components[0].size == 6
    assert s.components[0].hetero == ()
    assert s.components[0].aromatic is True
    assert s.fusions == ()


def test_naphthalene_fused():
    s = _ir("c1ccc2ccccc2c1")[0]
    assert len(s.components) == 2 and s.topology == "fused"
    assert len(s.atom_ids) == 10
    assert all(c.size == 6 and c.hetero == () for c in s.components)
    assert len(s.fusions) == 1
    assert len(s.fusions[0].shared) == 2
    assert s.fusions[0].bond is True


def test_indole_fused_hetero():
    s = _ir("c1ccc2[nH]ccc2c1")[0]
    assert len(s.components) == 2 and s.topology == "fused"
    assert len(s.atom_ids) == 9
    zs = {z for c in s.components for _, z in c.hetero}
    assert 7 in zs
    assert sorted(c.size for c in s.components) == [5, 6]
    assert len(s.fusions) == 1


def test_biphenyl_two_mono_systems():
    """Two unfused rings → two mono systems (not one fused)."""
    syss = _ir("c1ccc(-c2ccccc2)cc1")
    assert len(syss) == 2
    assert all(s.topology == "mono" for s in syss)
    assert all(len(s.components) == 1 for s in syss)


def test_open_chain_empty():
    assert _ir("CCCC") == []


def test_pyridine_hetero_mono():
    syss = _ir("c1ccncc1")
    assert len(syss) == 1
    s = syss[0]
    assert s.topology == "mono"
    assert any(z == 7 for c in s.components for _, z in c.hetero)


def test_fingerprint_placeholder():
    s = _ir("c1ccccc1")[0]
    assert s.fingerprint in ("", None)
