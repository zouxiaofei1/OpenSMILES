# IUPAC: P-66.5.1
# Layer: L1,L2,L4,L5
"""Simple acyclic mononitriles (alkanenitriles).

Nitrile carbon is C≡N; chain through that C. Retained acetonitrile;
C≥3 systematic …nitrile / …腈. Not alkyne (C≡C).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive
    ("C#N", "formonitrile", "甲腈"),
    ("CC#N", "acetonitrile", "乙腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("CCCCC#N", "pentanenitrile", "戊腈"),
    ("CCCCCC#N", "hexanenitrile", "己腈"),
    # negative
    ("C#CC", "propyne", "丙炔"),
    ("CCC", "propane", "丙烷"),
    ("CCN", "ethanamine", "乙胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_nitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
