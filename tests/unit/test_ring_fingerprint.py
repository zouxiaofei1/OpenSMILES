# IUPAC: P-25 / P-22
# Layer: L1
"""Ring layout fingerprint: sizes + fusion + hetero layout + aromatic."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_fingerprint import ring_fingerprint
from namepredict.layer1.ring_ir import RingSystemIR, build_ring_ir


def _ir(smiles: str) -> list[RingSystemIR]:
    mol = preprocess(smiles)
    assert mol is not None
    return build_ring_ir(mol)


def _fp(smiles: str) -> str:
    syss = _ir(smiles)
    assert len(syss) == 1
    fp = ring_fingerprint(syss[0])
    assert isinstance(fp, str) and fp
    return fp


def test_benzene_stable_mono_aromatic_size6():
    fp = _fp("c1ccccc1")
    assert "6" in fp
    assert "mono" in fp or "s=6" in fp
    assert fp == _fp("c1ccccc1")


def test_benzene_same_scaffold_same_fp():
    """Same layout scaffold always yields identical fingerprint."""
    assert _fp("c1ccccc1") == _fp("C1=CC=CC=C1")


def test_indole_ne_benzofuran():
    """N vs O hetero layout must differ (indole ≠ benzofuran)."""
    a = _fp("c1ccc2[nH]ccc2c1")
    b = _fp("c1ccc2occc2c1")
    assert a != b
    assert a and b


def test_indole_same_scaffold_same_fp():
    """Two representations of the same indole skeleton share fp."""
    assert _fp("c1ccc2[nH]ccc2c1") == _fp("c1cc2ccccc2[nH]1")


def test_naphthalene_fused_carbocycle():
    fp = _fp("c1ccc2ccccc2c1")
    assert fp != _fp("c1ccccc1")
    assert "6" in fp


def test_pyridine_ne_benzene():
    """Hetero mono layout differs from carbocycle mono."""
    assert _fp("c1ccncc1") != _fp("c1ccccc1")


def test_fingerprint_on_ir_when_wired():
    """If build_ring_ir fills fingerprint, it matches pure function."""
    s = _ir("c1ccccc1")[0]
    if s.fingerprint:
        assert s.fingerprint == ring_fingerprint(s)
