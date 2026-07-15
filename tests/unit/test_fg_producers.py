# IUPAC: P-44
# Layer: L2
"""FG parent producers registered in kind_registry (no hand-written _FG_TRY)."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.candidates import _fg_candidates

# Ordered producer __name__ snapshot (pre-migration candidates._FG_TRY).
_EXPECTED_FG_TRY_NAMES = (
    "_try_acid",
    "_try_anhydride",
    "_try_arene_other_fg",
    "_try_acyl_chloride",
    "_try_diester",
    "_try_ester",
    "_try_carbamate",
    "_try_carbonate",
    "_try_urea",
    "_try_guanidine",
    "_try_amide",
    "_try_nitrile",
    "_try_aldehyde",
    "_try_ketone",
    "_try_alcohol",
    "_try_thiol",
    "_try_benzenediamine",
    "_try_amine",
    "_try_hydrazine",
    "_try_phosphate",
    "_try_phosphonic",
    "_ether_parent",
    "_sulfide_parent",
    "_sulfoxide_parent",
    "_try_isocyanate",
    "_try_isothiocyanate",
    "_sulfonamide_parent",
    "_sulfonate_parent",
    "_sulfonyl_chloride_parent",
    "_sulfonic_acid_parent",
    "_try_boronic",
)

# Recent FG modules must be registered (by callable name substring).
_RECENT_FG_MARKERS = (
    "sulfonamide",
    "boronic",
    "urea",
    "guanidine",
    "hydrazine",
    "isocyanate",
    "sulfonate",
    "sulfonic_acid",
    "sulfoxide",
    "carbonate",
)

# Snapshot: FG kind sets before migration (must stay identical).
_KIND_CASES = [
    ("CC(=O)O", {"acid"}),
    ("CCO", {"alcohol"}),
    ("CS(=O)(=O)N", {"sulfonamide"}),
    ("OB(O)c1ccccc1", {"boronic"}),
    ("NC(=O)N", {"urea"}),
    ("NN", {"hydrazine"}),
    ("N=C(N)N", {"guanidine"}),
    ("CC(=O)Cl", {"acyl_chloride"}),
    ("CS(=O)(=O)O", {"sulfonic_acid"}),
]


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _kinds(cands: list[dict]) -> set[str]:
    return {c["kind"] for c in cands if c is not None}


def test_fg_try_fns_nonempty_and_order() -> None:
    fns = kr.fg_try_fns()
    assert len(fns) >= 20
    assert tuple(f.__name__ for f in fns) == _EXPECTED_FG_TRY_NAMES


def test_fg_try_fns_include_recent_fgs() -> None:
    names = " ".join(f.__name__ for f in kr.fg_try_fns())
    for marker in _RECENT_FG_MARKERS:
        assert marker in names, f"missing recent FG producer: {marker}"


def test_fg_try_fns_unique() -> None:
    names = [f.__name__ for f in kr.fg_try_fns()]
    assert len(names) == len(set(names))


def test_register_fg_try_appends() -> None:
    before = list(kr.fg_try_fns())

    def _try_dummy_fg(info: dict):
        return None

    kr.register_fg_try(_try_dummy_fg)
    assert kr.fg_try_fns()[-1] is _try_dummy_fg
    kr._FG_TRY.clear()
    for fn in before:
        kr.register_fg_try(fn)
    assert tuple(f.__name__ for f in kr.fg_try_fns()) == _EXPECTED_FG_TRY_NAMES


@pytest.mark.parametrize("smiles,expected", _KIND_CASES)
def test_fg_candidates_kinds_stable(smiles: str, expected: set[str]) -> None:
    assert _kinds(_fg_candidates(_info(smiles))) == expected


@pytest.mark.parametrize("smiles", ["C", "CCCC", "c1ccccc1"])
def test_fg_candidates_empty_or_no_fg_safe(smiles: str) -> None:
    """No FG / arene-only: must not raise; kinds may be empty."""
    cands = _fg_candidates(_info(smiles))
    assert isinstance(cands, list)
    assert all(isinstance(c, dict) and "kind" in c for c in cands)


def test_candidates_has_no_hand_fg_try_tuple() -> None:
    path = Path("src/namepredict/layer2/candidates.py")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id == "_FG_TRY":
                if isinstance(node.value, ast.Tuple) and len(node.value.elts) >= 5:
                    pytest.fail("hand-written _FG_TRY tuple still present")
