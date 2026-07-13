# IUPAC: P-22.1.1
# Layer: L2,L4,L5
"""Retained parent naphthalene (IUPAC P-22.1.1 / P-25): two fused aromatic
six-membered carbocycles. Unsubstituted + mono-methyl + mono-halo; lowest
locant (1- preferred over 2-). Fusion carbons are not substitution sites here.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2.naphthalene import _try_naphthalene_parent
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted (benchmark dual candidate)
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    # positive: mono-methyl — alpha (1) and beta (2)
    ("Cc1cccc2ccccc12", "1-methylnaphthalene", "1-甲基萘"),
    ("Cc1ccc2ccccc2c1", "2-methylnaphthalene", "2-甲基萘"),
    # positive: mono-halo
    ("Clc1cccc2ccccc12", "1-chloronaphthalene", "1-氯萘"),
    ("Clc1ccc2ccccc2c1", "2-chloronaphthalene", "2-氯萘"),
    # negative near-miss: must not become naphthalene
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CCCCCCCCCC", "decane", "癸烷"),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),  # indole retained parent (not naphthalene)
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_naphthalene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_decane() -> None:
    r = SMILESNNamer().name("c1ccc2ccccc2c1")
    assert r.success
    assert normalize_en(r.en) == "naphthalene"
    assert "decane" not in normalize_en(r.en)


def test_l2_parent_kind_unsub() -> None:
    mol = preprocess("c1ccc2ccccc2c1")
    assert mol is not None
    parent = _try_naphthalene_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "naphthalene"
    assert len(parent.get("chain") or []) == 10


def test_l2_rejects_benzene() -> None:
    mol = preprocess("c1ccccc1")
    assert mol is not None
    assert _try_naphthalene_parent(analyze(mol)) is None


def test_l2_rejects_indole() -> None:
    mol = preprocess("c1ccc2[nH]ccc2c1")
    assert mol is not None
    assert _try_naphthalene_parent(analyze(mol)) is None


def test_methyl_locant_alpha_is_1() -> None:
    r = SMILESNNamer().name("Cc1cccc2ccccc12")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-methylnaphthalene"
    assert "2-methyl" not in en
