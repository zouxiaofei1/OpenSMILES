# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained parent 1,3-benzothiazole (IUPAC P-22.2.1 / P-25 / P-62.2.1):
fused aromatic 6C+5(C3NS); S=1, N=3. Unsub / ≤2 halo+methyl+CF3;
2-amine with ≤1 halo or CF3 (gold: 6-(trifluoromethyl)-1,3-benzothiazol-2-amine).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.scaffold.benzazole import (
    _try_benzothiazole_parent,
    _try_benzothiazolamine_parent,
)
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted 1,3-benzothiazole
    ("c1nc2ccccc2s1", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("s1cnc2ccccc12", "1,3-benzothiazole", "1,3-苯并噻唑"),
    # positive: mono-halo / mono-methyl on ring
    ("Brc1ccc2ncsc2c1", "6-bromo-1,3-benzothiazole", "6-溴-1,3-苯并噻唑"),
    ("Cc1nc2ccccc2s1", "2-methyl-1,3-benzothiazole", "2-甲基-1,3-苯并噻唑"),
    # positive: unsubstituted 2-amine
    ("Nc1nc2ccccc2s1", "1,3-benzothiazol-2-amine", "2-氨基苯并噻唑"),
    # positive: dual gold — 6-(trifluoromethyl)-1,3-benzothiazol-2-amine
    (
        "FC(C1=CC2=C(N=C(S2)N)C=C1)(F)F",
        "6-(trifluoromethyl)-1,3-benzothiazol-2-amine",
        "2-氨基-6-三氟甲基苯并噻唑",
    ),
    # positive: 2-amine + mono-halo
    ("Nc1nc2ccc(Br)cc2s1", "6-bromo-1,3-benzothiazol-2-amine", "2-氨基-6-溴苯并噻唑"),
    # negative: near neighbors must not regress
    ("c1csc2ccccc12", "1-benzothiophene", "苯并[b]噻吩"),
    ("c1coc2ccccc12", "benzofuran", "苯并呋喃"),
    ("c1ccc2ncccc2c1", "quinoline", "喹啉"),
    ("c1cscn1", "1,3-thiazole", "噻唑"),
    ("c1ccccc1", "benzene", "苯"),
    ("Nc1ccccc1", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzothiazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_hexane() -> None:
    r = SMILESNNamer().name("c1nc2ccccc2s1")
    assert r.success
    assert normalize_en(r.en) == "1,3-benzothiazole"
    assert "hexane" not in normalize_en(r.en)


def test_gold_amine_not_methanamine() -> None:
    r = SMILESNNamer().name("FC(C1=CC2=C(N=C(S2)N)C=C1)(F)F")
    assert r.success
    en = normalize_en(r.en)
    assert en == "6-(trifluoromethyl)-1,3-benzothiazol-2-amine"
    assert "methanamine" not in en


def test_l2_parent_kind_unsub() -> None:
    mol = preprocess("c1nc2ccccc2s1")
    assert mol is not None
    parent = _try_benzothiazole_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzothiazole"
    assert len(parent.get("chain") or []) == 9
    assert parent.get("s_idx") is not None
    assert parent.get("n_idx") is not None


def test_l2_parent_kind_amine() -> None:
    mol = preprocess("Nc1nc2ccccc2s1")
    assert mol is not None
    parent = _try_benzothiazolamine_parent(analyze(mol))
    assert parent is not None
    assert parent.get("kind") == "benzothiazolamine"
    assert parent.get("amine_c_idx") is not None
    assert len(parent.get("chain") or []) == 9


def test_l2_rejects_benzothiophene() -> None:
    mol = preprocess("c1csc2ccccc12")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzothiazole_parent(info) is None
    assert _try_benzothiazolamine_parent(info) is None


def test_l2_rejects_thiazole() -> None:
    mol = preprocess("c1cscn1")
    assert mol is not None
    assert _try_benzothiazole_parent(analyze(mol)) is None


def test_l2_rejects_quinoline() -> None:
    mol = preprocess("c1ccc2ncccc2c1")
    assert mol is not None
    assert _try_benzothiazole_parent(analyze(mol)) is None


def test_simple_blocks_amine() -> None:
    """Primary 2-amine routes to benzothiazolamine, not plain benzothiazole."""
    mol = preprocess("Nc1nc2ccccc2s1")
    assert mol is not None
    info = analyze(mol)
    assert _try_benzothiazole_parent(info) is None
    assert _try_benzothiazolamine_parent(info) is not None


def test_bromo_locant_is_6() -> None:
    r = SMILESNNamer().name("Brc1ccc2ncsc2c1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "6-bromo-1,3-benzothiazole"
    assert "5-bromo" not in en


def test_amine_locant_is_2() -> None:
    r = SMILESNNamer().name("Nc1nc2ccccc2s1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1,3-benzothiazol-2-amine"
    assert "3-amine" not in en
