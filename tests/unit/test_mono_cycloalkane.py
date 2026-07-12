# IUPAC: P-22.1.1
# Layer: L1,L2,L5
"""Simple unsubstituted monocycloalkanes C3–C10 (P-22.1.1 / P-31).

cyclo + alkane stem; Chinese 环 + 烷. No side chains, no unsaturation.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted monocycloalkanes
    ("C1CC1", "cyclopropane", "环丙烷"),
    ("C1CCC1", "cyclobutane", "环丁烷"),
    ("C1CCCC1", "cyclopentane", "环戊烷"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("C1CCCCCC1", "cycloheptane", "环庚烷"),
    ("C1CCCCCCC1", "cyclooctane", "环辛烷"),
    ("C1CCCCCCCC1", "cyclononane", "环壬烷"),
    ("C1CCCCCCCCC1", "cyclodecane", "环癸烷"),
    # negative: acyclic / FG / unsaturated must not become cycloalkanes
    ("CCC", "propane", "丙烷"),
    ("CCCC", "butane", "丁烷"),
    ("C=C", "ethene", "乙烯"),
    ("C#C", "acetylene", "乙炔"),
    ("CCO", "ethanol", "乙醇"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCN", "ethanamine", "乙胺"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methylcyclohexane_not_cyclohexane() -> None:
    """Substituted ring is out of scope; must not be named cyclohexane."""
    r = SMILESNNamer().name("CC1CCCCC1")
    assert r.success
    assert normalize_en(r.en) != "cyclohexane"
