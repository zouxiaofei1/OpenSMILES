# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent 1-benzothiophene (IUPAC P-22.2.1 / P-25):
benzo[b]thiophene fused aromatic 6+5. S=1; unsubstituted, mono-methyl /
mono-halo; mono ring OH → 1-benzothiophen-n-ol / 苯并[b]噻吩-n-醇.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted 1-benzothiophene
    ("c1ccc2sccc2c1", "1-benzothiophene", "苯并[b]噻吩"),
    ("s1ccc2ccccc12", "1-benzothiophene", "苯并[b]噻吩"),
    # positive: dual gold — 5-fluoro-1-benzothiophene
    ("FC=1C=CC2=C(C=CS2)C1", "5-fluoro-1-benzothiophene", "5-氟苯并[b]噻吩"),
    # positive: dual gold — 1-benzothiophen-4-ol
    ("S1C=CC=2C1=CC=CC2O", "1-benzothiophen-4-ol", "苯并[b]噻吩-4-醇"),
    # positive: mono-methyl (α=2; benzo 5)
    ("Cc1cc2ccccc2s1", "2-methyl-1-benzothiophene", "2-甲基苯并[b]噻吩"),
    ("Cc1csc2ccccc12", "3-methyl-1-benzothiophene", "3-甲基苯并[b]噻吩"),
    # positive: mono-halo on benzo ring
    ("Clc1ccc2sccc2c1", "5-chloro-1-benzothiophene", "5-氯苯并[b]噻吩"),
    ("Brc1ccc2sccc2c1", "5-bromo-1-benzothiophene", "5-溴苯并[b]噻吩"),
    # negative: near neighbors must not regress
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("c1ccsc1", "thiophene", "噻吩"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzothiophene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2sccc2c1")
    assert r.success
    assert normalize_en(r.en) == "1-benzothiophene"
    assert "octane" not in normalize_en(r.en)


def test_ol_not_octanol() -> None:
    r = SMILESNNamer().name("S1C=CC=2C1=CC=CC2O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-benzothiophen-4-ol"
    assert "octan" not in en


def test_methyl_locant_is_2_not_3() -> None:
    """Cc1cc2ccccc2s1 is 2-methyl (α to S), not 3-methyl."""
    r = SMILESNNamer().name("Cc1cc2ccccc2s1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methyl-1-benzothiophene"
    assert "3-methyl" not in en


def test_fluoro_locant_is_5() -> None:
    """Benchmark dual: F on benzo ring is locant 5."""
    r = SMILESNNamer().name("FC=1C=CC2=C(C=CS2)C1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "5-fluoro-1-benzothiophene"
    assert "4-fluoro" not in en
