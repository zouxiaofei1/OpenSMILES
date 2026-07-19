# IUPAC: P-25 / P-63.1.4 / P-62.2.1
# Layer: L2,L4,L5
"""Generalized arene FG parent: mono/di OH/NH2 on fused aromatic rings.

Anchor molecules:
  Oc1ccc2ccccc2c1      → naphthalen-2-ol / 萘-2-酚
  Oc1cnc2c(O)cccc2c1   → quinoline-3,8-diol / 喹啉-3,8-二酚
  Nc1ccc2ccccc2c1      → naphthalen-2-amine / 萘-2-胺
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.arene_fg_parent import (
    _try_naphthalenediol_parent,
    _try_naphthalenol_parent,
    _try_naphthalenamine_parent,
    _try_quinolinediol_parent,
)
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # Anchor 1: mono OH on naphthalene
    ("Oc1ccc2ccccc2c1", "naphthalen-2-ol", "萘-2-酚"),
    # Anchor 3: mono NH2 on naphthalene
    ("Nc1ccc2ccccc2c1", "naphthalen-2-amine", "萘-2-胺"),
    # Anchor 2: di-OH on quinoline
    ("Oc1cnc2c(O)cccc2c1", "quinoline-3,8-diol", "喹啉-3,8-二酚"),
    # Negative: existing phenol/aniline must still work
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    # Negative: existing naphthalene must still work
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    # Negative: existing quinoline must still work
    ("c1ccc2ncccc2c1", "quinoline", "喹啉"),
    # Negative: existing quinolinol must still work
    ("Oc1ccc2ncccc2c1", "quinolin-6-ol", "喹啉-6-醇"),
    # Negative: existing pyridinol must still work
    ("Oc1ccccn1", "pyridin-2-ol", "吡啶-2-醇"),
    # Negative: chain alcohol must still work
    ("CCO", "ethanol", "乙醇"),
    # Negative: cyclohexanol must still work
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_fg_parent_full(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"Naming failed for {smiles}"
    assert normalize_en(r.en) == normalize_en(en), f"{smiles}: got {r.en!r}, expected {en!r}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh), f"{smiles}: zh got {r.zh!r}, expected {zh!r}"


def test_l2_naphthalenol_parent_kind() -> None:
    mol = preprocess("Oc1ccc2ccccc2c1")
    assert mol is not None
    parent = _try_naphthalenol_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "naphthalenol"
    assert len(parent.get("chain") or []) == 10
    assert parent.get("oh_c_idx") is not None


def test_l2_naphthalenamine_parent_kind() -> None:
    mol = preprocess("Nc1ccc2ccccc2c1")
    assert mol is not None
    parent = _try_naphthalenamine_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "naphthalenamine"
    assert len(parent.get("chain") or []) == 10
    assert parent.get("amine_c_idx") is not None


def test_l2_quinolinediol_parent_kind() -> None:
    mol = preprocess("Oc1cnc2c(O)cccc2c1")
    assert mol is not None
    parent = _try_quinolinediol_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "quinolinediol"
    assert len(parent.get("chain") or []) == 10
    assert parent.get("oh_c_idxs") is not None
    assert len(parent["oh_c_idxs"]) == 2


def test_l2_rejects_too_many_oh() -> None:
    """Tri-OH on naphthalene should not match (n_fg=2 cap currently)."""
    mol = preprocess("Oc1cc(O)c2c(O)cccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_naphthalenediol_parent(info) is None


def test_l2_rejects_chain_oh() -> None:
    """OH on side chain only, not on ring → should not match arene FG parent."""
    mol = preprocess("OCCc1ccc2ccccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_naphthalenol_parent(info) is None


def test_l2_rejects_bare_naphthalene() -> None:
    mol = preprocess("c1ccc2ccccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_naphthalenol_parent(info) is None
    assert _try_naphthalenamine_parent(info) is None
