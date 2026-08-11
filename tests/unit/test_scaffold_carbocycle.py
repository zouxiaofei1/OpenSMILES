# IUPAC: P-22.1 / P-31.1
# Layer: L2
"""Carbocycle ScaffoldSpec registration: cycloalkane / ene / polyene."""
from __future__ import annotations

import pytest

from namepredict.layer2.scaffold.specs import CARBOCYCLE_SPECS, ScaffoldSpec, get_spec


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
