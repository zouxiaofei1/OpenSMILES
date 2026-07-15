# IUPAC: P-31.1 / P-44
# Layer: L2
"""Unsat parent producers registered in kind_registry (no hand-written unsat tries)."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.candidates import _unsat_candidates

# Ordered producer __name__ snapshot (pre-migration candidates._unsat_candidates).
_EXPECTED_UNSAT_TRY_NAMES = (
    "_try_alkyne",
    "_try_polyene",
    "_try_alkene",
)

# Snapshot: unsat kind sets before migration (must stay identical).
_KIND_CASES = [
    ("C=C", {"alkene"}),
    ("C#C", {"alkyne"}),
    ("C=CC=C", {"polyene"}),
    ("CC=C", {"alkene"}),
    ("CC#C", {"alkyne"}),
    ("C=CC=CC=C", {"polyene"}),
]


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _kinds(cands: list[dict]) -> set[str]:
    return {c["kind"] for c in cands if c is not None}


def test_unsat_try_fns_order() -> None:
    fns = kr.unsat_try_fns()
    assert tuple(f.__name__ for f in fns) == _EXPECTED_UNSAT_TRY_NAMES


def test_unsat_try_fns_unique() -> None:
    names = [f.__name__ for f in kr.unsat_try_fns()]
    assert len(names) == len(set(names))


def test_register_unsat_try_appends() -> None:
    before = list(kr.unsat_try_fns())

    def _try_dummy_unsat(info: dict):
        return None

    kr.register_unsat_try(_try_dummy_unsat)
    assert kr.unsat_try_fns()[-1] is _try_dummy_unsat
    kr._UNSAT_TRY.clear()
    for fn in before:
        kr.register_unsat_try(fn)
    assert tuple(f.__name__ for f in kr.unsat_try_fns()) == _EXPECTED_UNSAT_TRY_NAMES


@pytest.mark.parametrize("smiles,expected", _KIND_CASES)
def test_unsat_candidates_kinds_stable(smiles: str, expected: set[str]) -> None:
    assert _kinds(_unsat_candidates(_info(smiles))) == expected


@pytest.mark.parametrize("smiles", ["CCO", "C", "CCCC", "CC(=O)O"])
def test_unsat_candidates_empty_for_non_unsat(smiles: str) -> None:
    """Alcohol / alkane / acid: unsat-only path must not yield parents."""
    cands = _unsat_candidates(_info(smiles))
    assert cands == []


def test_unsat_candidates_match_registry_iteration() -> None:
    """candidates._unsat_candidates kinds match direct unsat_try_fns() iteration."""
    for smiles in ("C=C", "C#C", "C=CC=C", "CCO"):
        info = _info(smiles)
        via_cands = _kinds(_unsat_candidates(info))
        via_reg = _kinds(
            [c for fn in kr.unsat_try_fns() if (c := fn(info)) is not None]
        )
        assert via_cands == via_reg


def test_candidates_has_no_hand_unsat_try_tuple() -> None:
    path = Path("src/namepredict/layer2/candidates.py")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id in ("_UNSAT_TRY", "_UNSAT_PRODUCERS"):
                if isinstance(node.value, ast.Tuple) and len(node.value.elts) >= 2:
                    pytest.fail("hand-written unsat try tuple still in candidates")
