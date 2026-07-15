# IUPAC: P-22.1 / P-31.1
# Layer: L2
"""Carbocycle ScaffoldSpec + builder: cycloalkane / ene / polyene recognition."""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2.scaffold.builders.carbocycle import try_carbocycle
from namepredict.layer2.scaffold.specs import (
    CARBOCYCLE_SPECS,
    ScaffoldSpec,
    get_spec,
)

# Positive: mono carbocycle recognition (not namer).
CASES = [
    ("C1CCCCC1", "cycloalkane", 6, 0),
    ("C1=CCCCC1", "cycloalkene", 6, 1),
    ("C1=CC=CCC1", "cyclopolyene", 6, 2),  # 1,3
    ("C1C=CC=CC1", "cyclopolyene", 6, 2),  # 1,3 alt
    ("C1=CCC=CC1", "cyclopolyene", 6, 2),  # 1,4
    ("C1CCC=CC1", "cycloalkene", 6, 1),
    ("C1CCCC1", "cycloalkane", 5, 0),
    ("C1=CCCC1", "cycloalkene", 5, 1),
]

# Negative: aromatic arene / open polyene — not carbocycle scaffold.
NEG = [
    "c1ccccc1",  # benzene (arene owns)
    "C=CC=C",  # open chain polyene
    "c1ccncc1",  # hetero aromatic
]


def _mol(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    return mol


@pytest.mark.parametrize("smiles,spec_id,size,n_double", CASES)
def test_try_carbocycle_hit(smiles, spec_id, size, n_double):
    hit = try_carbocycle(_mol(smiles))
    assert hit is not None
    assert hit.spec_id == spec_id
    assert hit.ring_size == size
    assert hit.n_double == n_double
    assert len(hit.atom_ids) == size
    assert get_spec(hit.spec_id) is not None


@pytest.mark.parametrize("smiles", NEG)
def test_try_carbocycle_reject(smiles):
    assert try_carbocycle(_mol(smiles)) is None


def test_carbocycle_specs_registered():
    for sid in ("cycloalkane", "cycloalkene", "cyclopolyene"):
        sp = get_spec(sid)
        assert isinstance(sp, ScaffoldSpec)
        assert sp.id == sid
        assert sp.ring == "carbo"
        assert sp.n_rings == 1
        assert sp in CARBOCYCLE_SPECS or sp.id in {s.id for s in CARBOCYCLE_SPECS}


def test_numbering_modes():
    assert get_spec("cycloalkane").numbering.mode == "carbocycle_free"
    assert get_spec("cycloalkene").numbering.mode == "carbocycle_free"
    assert get_spec("cyclopolyene").numbering.mode == "poly_unsat"
