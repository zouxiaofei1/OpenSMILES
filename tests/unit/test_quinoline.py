# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent quinoline / isoquinoline (IUPAC P-22.2.1 / P-25):
fused aromatic 6+6, N=1 (quinoline) or N=2 (isoquinoline).
Unsub, ≤2 halo/methyl; mono ring OH → quinolin-n-ol; mono COOH.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.quinoline import (
    _try_isoquinoline_parent,
    _try_quinoline_parent,
    _try_quinolinecarboxylic_parent,
    _try_quinolinol_parent,
)
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted
    ("c1ccc2ncccc2c1", "quinoline", "喹啉"),
    ("c1ccc2cnccc2c1", "isoquinoline", "异喹啉"),
    # positive: dual gold — 2-chloroquinolin-6-ol
    ("ClC1=NC2=CC=C(C=C2C=C1)O", "2-chloroquinolin-6-ol", "2-氯喹啉-6-醇"),
    # positive: dual gold — 2-chloro-4-methylquinolin-6-ol
    ("ClC1=NC2=CC=C(C=C2C(=C1)C)O", "2-chloro-4-methylquinolin-6-ol", "2-氯-4-甲基喹啉-6-醇"),
    # positive: mono COOH (en gold)
    ("O=C(O)c1ccc2ccccc2n1", "quinoline-2-carboxylic acid", "喹啉-2-甲酸"),
    # positive: mono-halo / mono-methyl
    ("Clc1ccc2ncccc2c1", "6-chloroquinoline", "6-氯喹啉"),
    ("Cc1ccc2ncccc2c1", "6-methylquinoline", "6-甲基喹啉"),
    ("Clc1cncc2ccccc12", "4-chloroisoquinoline", "4-氯异喹啉"),
    # negative: near neighbors must not regress
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1ccc2sccc2c1", "1-benzothiophene", "苯并[b]噻吩"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_quinoline_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_decane() -> None:
    r = SMILESNNamer().name("c1ccc2ncccc2c1")
    assert r.success
    assert normalize_en(r.en) == "quinoline"
    assert "decane" not in normalize_en(r.en)
    assert "naphthalene" not in normalize_en(r.en)


def test_iso_not_quinoline() -> None:
    r = SMILESNNamer().name("c1ccc2cnccc2c1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "isoquinoline"
    assert en != "quinoline"


def test_ol_not_decanol() -> None:
    r = SMILESNNamer().name("ClC1=NC2=CC=C(C=C2C=C1)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-chloroquinolin-6-ol"
    assert "decan" not in en
    assert "hexan" not in en


def test_l2_parent_kind_quinoline() -> None:
    mol = preprocess("c1ccc2ncccc2c1")
    assert mol is not None
    parent = _try_quinoline_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "quinoline"
    assert len(parent.get("chain") or []) == 10
    assert parent.get("n_idx") is not None


def test_l2_parent_kind_isoquinoline() -> None:
    mol = preprocess("c1ccc2cnccc2c1")
    assert mol is not None
    parent = _try_isoquinoline_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "isoquinoline"
    assert len(parent.get("chain") or []) == 10
    assert parent.get("n_idx") is not None


def test_l2_parent_kind_ol() -> None:
    mol = preprocess("ClC1=NC2=CC=C(C=C2C=C1)O")
    assert mol is not None
    parent = _try_quinolinol_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "quinolinol"
    assert parent.get("oh_c_idx") is not None
    assert len(parent.get("chain") or []) == 10


def test_l2_parent_kind_carboxylic() -> None:
    mol = preprocess("O=C(O)c1ccc2ccccc2n1")
    assert mol is not None
    parent = _try_quinolinecarboxylic_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "quinolinecarboxylic"
    assert parent.get("cooh_c_idx") is not None
    assert parent.get("ring_attach_idx") is not None


def test_l2_rejects_naphthalene() -> None:
    mol = preprocess("c1ccc2ccccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is None
    assert _try_isoquinoline_parent(info) is None


def test_l2_rejects_pyridine() -> None:
    mol = preprocess("c1ccncc1")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is None
    assert _try_isoquinoline_parent(info) is None


def test_l2_rejects_indole() -> None:
    mol = preprocess("c1ccc2[nH]ccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is None
    assert _try_isoquinoline_parent(info) is None


def test_l2_rejects_benzofuran() -> None:
    mol = preprocess("c1ccc2occc2c1")
    assert mol is not None
    assert _try_quinoline_parent(analyze(mol)) is None


def test_simple_blocks_ol() -> None:
    """Ring OH routes to quinolinol, not plain quinoline."""
    mol = preprocess("ClC1=NC2=CC=C(C=C2C=C1)O")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is None
    assert _try_quinolinol_parent(info) is not None


def test_simple_blocks_carboxylic() -> None:
    mol = preprocess("O=C(O)c1ccc2ccccc2n1")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is None
    assert _try_quinolinecarboxylic_parent(info) is not None


def test_quinoline_not_iso() -> None:
    mol = preprocess("c1ccc2ncccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_quinoline_parent(info) is not None
    assert _try_isoquinoline_parent(info) is None


def test_iso_not_plain_quinoline_parent() -> None:
    mol = preprocess("c1ccc2cnccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_isoquinoline_parent(info) is not None
    assert _try_quinoline_parent(info) is None


def test_ol_locants_2_and_6() -> None:
    r = SMILESNNamer().name("ClC1=NC2=CC=C(C=C2C=C1)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-chloroquinolin-6-ol"
    assert "quinolin-5-ol" not in en
    assert "quinolin-7-ol" not in en


def test_disub_ol_locants() -> None:
    r = SMILESNNamer().name("ClC1=NC2=CC=C(C=C2C(=C1)C)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-chloro-4-methylquinolin-6-ol"
    assert "3-methyl" not in en
