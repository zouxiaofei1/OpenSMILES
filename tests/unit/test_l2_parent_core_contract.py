# IUPAC: P-44 / architecture L2 parent hub
# Layer: L2
"""L2 parent_core is the sole authority for parent assembly / chain / gate helpers.

Producer modules must not import private helpers from parent_selector.
parent_selector keeps FG try + select_parent; helpers live in parent_core.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_ROOT = Path(__file__).resolve().parents[2]
_L2 = _ROOT / "src" / "namepredict" / "layer2"

# Helpers that must be sourced from parent_core, not parent_selector.
_HELPERS = frozenset({
    "_parent_dict",
    "_no_fgs",
    "_longest_from",
    "_longest_chain",
    "_c_idxs",
})

# Modules allowed to re-export / host try surface.
_EXEMPT = frozenset({
    "parent_selector.py",
    "parent_core.py",
    "fg_helpers.py",
    "__init__.py",
})


def _parent_selector_helper_hits(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom) or not node.module:
            continue
        if node.module != "namepredict.layer2.parent_selector":
            continue
        for alias in node.names:
            if alias.name in _HELPERS:
                hits.append(f"{alias.name} (L{node.lineno})")
    return hits


def _producer_paths() -> list[Path]:
    return sorted(
        p for p in _L2.rglob("*.py")
        if p.name not in _EXEMPT and p.is_file()
    )


def test_parent_core_module_exists() -> None:
    assert (_L2 / "parent_core.py").is_file()


def test_parent_core_exports_helpers() -> None:
    import namepredict.layer2.parent_core as core

    missing = [h for h in sorted(_HELPERS) if not hasattr(core, h)]
    assert missing == [], f"parent_core missing helpers: {missing}"


@pytest.mark.parametrize("path", _producer_paths(), ids=lambda p: p.relative_to(_L2).as_posix())
def test_producers_do_not_import_helpers_from_parent_selector(path: Path) -> None:
    hits = _parent_selector_helper_hits(path)
    assert hits == [], (
        f"{path.relative_to(_L2).as_posix()} must not import helpers "
        f"from parent_selector: {hits}"
    )


# End-to-end bilingual smoke (zero-behavior after cut).
_E2E = [
    # acid / chain
    ("CC(=O)O", "acetic acid", "乙酸"),
    # ester
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_parent_core_e2e_bilingual(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_alkyl_acid_negative_not_benzene() -> None:
    """Negative: simple aliphatic acid stays alkyl (not benzoic)."""
    r = SMILESNNamer().name("CCC(=O)O")
    assert r.success
    assert normalize_en(r.en) == "propanoic acid"
    assert "benzoic" not in normalize_en(r.en)
    assert "苯" not in normalize_zh(r.zh)
