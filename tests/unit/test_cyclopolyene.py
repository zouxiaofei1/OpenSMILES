# IUPAC: P-31.1 / P-22.1.3
# Layer: L2,L4,L5
"""Unsubstituted monocyclic polyenes (cyclopolyene): cyclohexadiene etc.

Parent = mono all-C carbocycle with ≥2 endocyclic non-aromatic C=C; no main FG.
Ene locant set lowest (P-31.1); EN cyclohexa-1,3-diene; ZH 环己-1,3-二烯.
Mono cycloalkene and open polyene must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: ring dienes
    ("C1=CC=CCC1", "cyclohexa-1,3-diene", "环己-1,3-二烯"),
    ("C1C=CC=CC1", "cyclohexa-1,3-diene", "环己-1,3-二烯"),
    ("C1=CCC=CC1", "cyclohexa-1,4-diene", "环己-1,4-二烯"),
    ("C1=CC=CC1", "cyclopenta-1,3-diene", "环戊-1,3-二烯"),
    # negative: mono cycloalkene / open polyene must not become cyclopolyene wrong
    ("C1=CCCCC1", "cyclohexene", "环己烯"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_cyclopolyene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C1=CC=CCC1", "C1=CCC=CC1", "C1=CC=CC1"])
def test_not_methane_or_alkane(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert en != "methane"
    assert "diene" in en
    assert "二烯" in (r.zh or "")
