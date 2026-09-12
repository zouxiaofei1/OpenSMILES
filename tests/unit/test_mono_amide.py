# IUPAC: P-66.1.1
# Layer: L1,L2,L4,L5
"""Simple primary unsubstituted alkanamides (–CONH2).

Primary amide: carbonyl C with =O and N whose only carbon neighbor is that C.
Retained formamide/acetamide; C≥3 systematic …amide / …酰胺.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: retained + straight-chain mono primary amides
    ("C(=O)N", "formamide", "甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCC(=O)N", "propanamide", "丙酰胺"),
    ("CCCC(=O)N", "butanamide", "丁酰胺"),
    ("CCCCC(=O)N", "pentanamide", "戊酰胺"),
    ("CCCCCC(=O)N", "hexanamide", "己酰胺"),
    # alkene / cycloalkane must not become amides
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("CCN", "ethanamine", "乙胺"),
    ("CCC", "propane", "丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
