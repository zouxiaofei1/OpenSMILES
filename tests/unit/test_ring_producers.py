# IUPAC: P-44 / architecture
# Layer: L2
"""Ring parent producers registered in kind_registry (no hand-written _RING_TRY)."""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.candidates import _ring_candidates
from namepredict.namer import SMILESNNamer

# Expected ordered producer __name__ snapshot (ring_producers._RING_PRODUCERS).
_EXPECTED_RING_TRY_NAMES = (
    "_try_anthraquinone_parent",
    "_try_anthracene_parent",
    "_try_quinazoline_parent",
    "_try_quinoxaline_parent",
    "_try_naphthalene_parent",
    "_try_indole_parent",
    "_try_indazole_parent",
    "_try_benzofuran_parent",
    "_try_benzothiophene_parent",
    "_try_benzothiazole_parent",
    "_try_benzoxazole_parent",
    "_try_benzimidazolamine_parent",
    "_try_benzimidazole_parent",
    "_try_quinoline_parent",
    "_try_isoquinoline_parent",
    "_try_pyridine_parent",
    "_try_diazine_parent",
    "_try_imidazole_parent",
    "_try_pyrazole_parent",
    "_try_azole13_parent",
    "_try_hetero5_parent",
    "_try_sat_hetero_parent",
    "try_sat_hetero_repl",
    "_try_simple_benzene",
    "_try_simple_cycloalkane",
    "_try_simple_cyclopolyene",
)


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def test_ring_try_fns_count_and_order() -> None:
    fns = kr.ring_try_fns()
    assert len(fns) == len(_EXPECTED_RING_TRY_NAMES)
    assert tuple(f.__name__ for f in fns) == _EXPECTED_RING_TRY_NAMES


def test_ring_try_fns_unique() -> None:
    names = [f.__name__ for f in kr.ring_try_fns()]
    assert len(names) == len(set(names))


def test_register_ring_try_appends() -> None:
    """register_ring_try is public; test via temporary fn then restore."""
    before = list(kr.ring_try_fns())

    def _try_dummy_parent(info: dict):
        return None

    kr.register_ring_try(_try_dummy_parent)
    assert kr.ring_try_fns()[-1] is _try_dummy_parent
    # restore (registry list is mutable by design for bootstrap)
    kr._RING_TRY.clear()
    for fn in before:
        kr.register_ring_try(fn)
    assert tuple(f.__name__ for f in kr.ring_try_fns()) == _EXPECTED_RING_TRY_NAMES


@pytest.mark.parametrize(
    "smiles,kind",
    [
        ("c1ccccc1", "benzene"),
        ("c1ccncc1", "pyridine"),
        ("C1CCOC1", "oxolane"),
        ("c1ccc2ccccc2c1", "naphthalene"),
        ("c1ccc2cc3ccccc3cc2c1", "anthracene"),
    ],
)
def test_ring_candidates_via_registry(smiles: str, kind: str) -> None:
    kinds = {c["kind"] for c in _ring_candidates(_info(smiles))}
    assert kind in kinds


def test_candidates_has_no_hand_ring_try_tuple() -> None:
    """candidates.py must not keep a large hand-written _RING_TRY name tuple."""
    path = Path("src/namepredict/layer2/candidates.py")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if isinstance(t, ast.Name) and t.id == "_RING_TRY":
                # allow alias to registry call, not a big Tuple of Names
                if isinstance(node.value, ast.Tuple) and len(node.value.elts) >= 5:
                    pytest.fail("hand-written _RING_TRY tuple still present")


# end-to-end retained parents still assemble
_E2E = [
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_e2e_ring_parents(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unknown_not_a_producer_name() -> None:
    names = {f.__name__ for f in kr.ring_try_fns()}
    assert "_try_no_such_parent" not in names
