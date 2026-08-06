# IUPAC: P-65 / architecture L4 orient table
# Layer: L4
"""L4 cycloalkane* kinds must use an explicit orienter table (no startswith fallback).

Known Round D exo FG kinds map to _orient_benzoic (ring_attach_idx).
Unknown cycloalkane* kinds must fall through to default alkane orient —
not silently reuse the benzoic ring-fixed path.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.layer4 import numbering as num

_ROOT = Path(__file__).resolve().parents[2]
_NUMBERING = _ROOT / "src" / "namepredict" / "layer4" / "numbering.py"

# Known cycloalkane* kinds that must be explicit in _kind_orienters().
# Positive cases: table membership + expected orienter function.
_KNOWN_CYCLOALKANE_KINDS = (
    # base / poly FG (already explicit before this round)
    "cycloalkane",
    "cycloalkanediol",
    "cycloalkanedione",
    "cycloalkane_polycarboxylic",
    # Only cycloalkanecarboxylic is wired (typed_acid_kind derivation);
    # carbaldehyde/carbonitrile/etc. exo FG kinds are un-wired and not produced.
    "cycloalkanecarboxylic",
)

# Exo FG kinds that must share _orient_benzoic (ring_attach_idx fixed).
_EXO_FG_KINDS = (
    "cycloalkanecarboxylic",
)


def test_no_cycloalkane_startswith_fallback() -> None:
    """startswith('cycloalkane') orient fallback must be deleted."""
    tree = ast.parse(_NUMBERING.read_text(encoding="utf-8"), filename=str(_NUMBERING))
    hits: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not isinstance(func, ast.Attribute) or func.attr != "startswith":
            continue
        if not node.args:
            continue
        arg0 = node.args[0]
        if isinstance(arg0, ast.Constant) and str(arg0.value).startswith("cycloalkane"):
            hits.append(f"L{node.lineno}")
    assert hits == [], f"startswith cycloalkane fallback still present: {hits}"


@pytest.mark.parametrize("kind", _KNOWN_CYCLOALKANE_KINDS)
def test_known_cycloalkane_kind_in_explicit_table(kind: str) -> None:
    table = num._kind_orienters()
    assert kind in table, f"{kind!r} missing from explicit orienter table"


@pytest.mark.parametrize("kind", _EXO_FG_KINDS)
def test_exo_fg_kinds_use_orient_benzoic(kind: str) -> None:
    table = num._kind_orienters()
    assert table.get(kind) is num._orient_benzoic, (
        f"{kind!r} must map to _orient_benzoic, got {table.get(kind)}"
    )


def test_unknown_cycloalkane_kind_not_benzoic() -> None:
    """Negative: unknown cycloalkane_foo must not silently use benzoic orient.

    With ring_attach_idx present, benzoic rotates attach to front; alkane default
    does not. After removing startswith, unknown kind → _orient_alkane.
    """
    chain = [10, 11, 12, 13]
    parent = {"kind": "cycloalkane_foo", "ring_attach_idx": 12}
    subs: list = []
    out = num._orient_by_kind("cycloalkane_foo", chain, parent, subs)
    # Default alkane with empty subs returns chain unchanged.
    assert out == chain
    # Benzoic would rotate 12 to front → [12, 13, 10, 11].
    assert out != num._orient_benzoic(chain, parent, subs)


def test_known_exo_still_uses_ring_attach() -> None:
    """Regression: known exo kind still fixes ring_attach_idx as locant 1."""
    chain = [10, 11, 12, 13]
    parent = {"kind": "cycloalkanecarboxylic", "ring_attach_idx": 12}
    out = num._orient_by_kind("cycloalkanecarboxylic", chain, parent, [])
    assert out[0] == 12
