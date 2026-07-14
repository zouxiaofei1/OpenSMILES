# IUPAC: P-22.2.1 / architecture
# Layer: L2
"""Data-driven fused 5+6 mono-hetero engine (benzofuran / benzothiophene)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.fused56 import (
    Fused56MonoSpec,
    _mono_parts,
    _try_mono_fused56,
)
from namepredict.namer import SMILESNNamer

BF = Fused56MonoSpec(
    kind="benzofuran", hetero_z=8, hetero_key="o_idx", sub_cap=1,
)
BT = Fused56MonoSpec(
    kind="benzothiophene", hetero_z=16, hetero_key="s_idx", sub_cap=1,
)


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def test_bf_parts_unsubstituted() -> None:
    p = _mono_parts(_info("c1ccc2occc2c1"), BF)
    assert p is not None
    five, six, h, ba, bb = p
    assert len(five) == 5 and len(six) == 6
    assert h is not None


def test_bt_parts_unsubstituted() -> None:
    p = _mono_parts(_info("c1ccc2sccc2c1"), BT)
    assert p is not None


def test_bf_try_parent() -> None:
    parent = _try_mono_fused56(_info("c1ccc2occc2c1"), BF)
    assert parent is not None
    assert parent["kind"] == "benzofuran"
    assert "o_idx" in parent
    assert len(parent["chain"]) == 9


def test_bt_try_parent() -> None:
    parent = _try_mono_fused56(_info("c1ccc2sccc2c1"), BT)
    assert parent is not None
    assert parent["kind"] == "benzothiophene"
    assert "s_idx" in parent


def test_furan_not_benzofuran() -> None:
    assert _try_mono_fused56(_info("c1ccoc1"), BF) is None


def test_indole_not_benzofuran() -> None:
    assert _try_mono_fused56(_info("c1ccc2[nH]ccc2c1"), BF) is None


# e2e: public try_* still work after thin wrappers
_E2E = [
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1ccc2sccc2c1", "1-benzothiophene", "苯并[b]噻吩"),
    ("Cc1cc2ccccc2o1", "2-methylbenzofuran", "2-甲基苯并呋喃"),
    ("FC=1C=CC2=C(C=CS2)C1", "5-fluoro-1-benzothiophene", "5-氟苯并[b]噻吩"),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_e2e_fused56_mono(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
