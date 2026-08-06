# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent benzofuran (IUPAC P-22.2.1 / P-25):
benzo[b]furan fused aromatic 6+5. O=1; unsubstituted, mono-methyl /
mono-halo; mono primary amine → benzofuran-n-amine / 苯并呋喃-n-胺.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.fused56_mono import (
    _try_benzofuran_parent,
    _try_benzofuranamine_parent,
)
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted benzofuran
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("o1ccc2ccccc12", "benzofuran", "苯并呋喃"),
    # positive: dual gold — benzofuran-2-amine
    ("O1C(=CC2=C1C=CC=C2)N", "benzofuran-2-amine", "苯并呋喃-2-胺"),
    # positive: mono-methyl (α=2, β=3; benzo 5)
    ("Cc1cc2ccccc2o1", "2-methylbenzofuran", "2-甲基苯并呋喃"),
    ("Cc1coc2ccccc12", "3-methylbenzofuran", "3-甲基苯并呋喃"),
    # positive: mono-halo on benzo ring
    ("Clc1ccc2occc2c1", "5-chlorobenzofuran", "5-氯苯并呋喃"),
    ("Brc1ccc2occc2c1", "5-bromobenzofuran", "5-溴苯并呋喃"),
    # negative: near neighbors must not regress
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("c1ccoc1", "furan", "呋喃"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzofuran_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2occc2c1")
    assert r.success
    assert normalize_en(r.en) == "benzofuran"
    assert "octane" not in normalize_en(r.en)


def test_amine_not_octanamine() -> None:
    r = SMILESNNamer().name("O1C(=CC2=C1C=CC=C2)N")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzofuran-2-amine"
    assert "octan" not in en


def test_l2_parent_kind_unsub() -> None:
    mol = preprocess("c1ccc2occc2c1")
    assert mol is not None
    parent = _try_benzofuran_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzofuran"
    assert len(parent.get("chain") or []) == 9
    assert parent.get("o_idx") is not None


def test_l2_parent_kind_amine() -> None:
    mol = preprocess("O1C(=CC2=C1C=CC=C2)N")
    assert mol is not None
    parent = _try_benzofuranamine_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzofuranamine"
    assert parent.get("amine_c_idx") is not None
    assert len(parent.get("chain") or []) == 9


def test_l2_rejects_indole() -> None:
    mol = preprocess("c1ccc2[nH]ccc2c1")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzofuran_parent(info) is None
    assert _try_benzofuranamine_parent(info) is None


def test_l2_rejects_furan() -> None:
    mol = preprocess("c1ccoc1")
    assert mol is not None
    assert _try_benzofuran_parent(analyze(mol)) is None


def test_l2_rejects_naphthalene() -> None:
    mol = preprocess("c1ccc2ccccc2c1")
    assert mol is not None
    assert _try_benzofuran_parent(analyze(mol)) is None


def test_simple_blocks_amine() -> None:
    """Primary amine routes to benzofuranamine, not plain benzofuran."""
    mol = preprocess("O1C(=CC2=C1C=CC=C2)N")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzofuran_parent(info) is None
    assert _try_benzofuranamine_parent(info) is not None


def test_methyl_locant_is_2_not_3() -> None:
    """Cc1cc2ccccc2o1 is 2-methyl (α to O), not 3-methyl."""
    r = SMILESNNamer().name("Cc1cc2ccccc2o1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylbenzofuran"
    assert "3-methyl" not in en
