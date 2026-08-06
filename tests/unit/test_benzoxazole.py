# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent 1,3-benzoxazole (IUPAC P-22.2.1 / P-25 / P-62.2.1):
fused aromatic 6C+5(C3NO); O=1, N=3 (1,3-oxazole type). Unsub / ≤2
halo+methyl+CF3; C2 primary amine → 1,3-benzoxazol-2-amine.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.benzoxazole import (
    _try_benzoxazole_parent,
    _try_benzoxazolamine_parent,
)
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted 1,3-benzoxazole
    ("c1ccc2ocnc2c1", "1,3-benzoxazole", "1,3-苯并噁唑"),
    ("o1cnc2ccccc12", "1,3-benzoxazole", "1,3-苯并噁唑"),
    # positive: mono-halo / mono-methyl on ring
    ("Brc1ccc2ncoc2c1", "6-bromo-1,3-benzoxazole", "6-溴-1,3-苯并噁唑"),
    ("Cc1nc2ccccc2o1", "2-methyl-1,3-benzoxazole", "2-甲基-1,3-苯并噁唑"),
    # positive: unsubstituted 2-amine
    ("Nc1nc2ccccc2o1", "1,3-benzoxazol-2-amine", "2-氨基苯并噁唑"),
    # positive: 2-amine + mono-halo
    ("Nc1nc2ccc(Br)cc2o1", "6-bromo-1,3-benzoxazol-2-amine", "2-氨基-6-溴苯并噁唑"),
    # negative: near neighbors must not regress
    ("c1nc2ccccc2s1", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("c1ccc2[nH]cnc2c1", "1H-benzimidazole", "1H-苯并咪唑"),
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("c1coc2ccccc12", "benzofuran", "苯并呋喃"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzoxazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_hexane() -> None:
    r = SMILESNNamer().name("c1ccc2ocnc2c1")
    assert r.success
    assert normalize_en(r.en) == "1,3-benzoxazole"
    assert "hexane" not in normalize_en(r.en)


def test_amine_not_methanamine() -> None:
    r = SMILESNNamer().name("Nc1nc2ccccc2o1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1,3-benzoxazol-2-amine"
    assert "methanamine" not in en


def test_l2_parent_kind_unsub() -> None:
    mol = preprocess("c1ccc2ocnc2c1")
    assert mol is not None
    parent = _try_benzoxazole_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzoxazole"
    assert len(parent.get("chain") or []) == 9
    assert parent.get("o_idx") is not None
    assert parent.get("n_idx") is not None


def test_l2_parent_kind_amine() -> None:
    mol = preprocess("Nc1nc2ccccc2o1")
    assert mol is not None
    parent = _try_benzoxazolamine_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzoxazolamine"
    assert parent.get("amine_c_idx") is not None
    assert len(parent.get("chain") or []) == 9


def test_l2_rejects_benzothiazole() -> None:
    mol = preprocess("c1nc2ccccc2s1")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzoxazole_parent(info) is None
    assert _try_benzoxazolamine_parent(info) is None


def test_l2_rejects_benzimidazole() -> None:
    mol = preprocess("c1ccc2[nH]cnc2c1")
    assert mol is not None
    assert _try_benzoxazole_parent(analyze(mol)) is None


def test_l2_rejects_anisole() -> None:
    mol = preprocess("COc1ccccc1")
    assert mol is not None
    assert _try_benzoxazole_parent(analyze(mol)) is None


def test_simple_blocks_amine() -> None:
    """Primary 2-amine routes to benzoxazolamine, not plain benzoxazole."""
    mol = preprocess("Nc1nc2ccccc2o1")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzoxazole_parent(info) is None
    assert _try_benzoxazolamine_parent(info) is not None


def test_bromo_locant_is_6() -> None:
    r = SMILESNNamer().name("Brc1ccc2ncoc2c1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "6-bromo-1,3-benzoxazole"
    assert "5-bromo" not in en


def test_amine_locant_is_2() -> None:
    r = SMILESNNamer().name("Nc1nc2ccccc2o1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1,3-benzoxazol-2-amine"
    assert "3-amine" not in en
