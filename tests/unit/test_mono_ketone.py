# IUPAC: P-64.2.1
# Layer: L1,L2,L4,L5
"""Simple acyclic monoketones (alkanones): stem + locant + one/酮.

Carbonyl carbon has exactly two carbon neighbors and is not carboxyl.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: straight-chain monoketones
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("CCC(CC)=O", "pentan-3-one", "戊-3-酮"),
    ("CCCC(=O)CC", "hexan-3-one", "己-3-酮"),
    ("CCCC(=O)C", "pentan-2-one", "戊-2-酮"),
    ("CCCCC(=O)CC", "heptan-3-one", "庚-3-酮"),
    ("CCCCCC(=O)C", "heptan-2-one", "庚-2-酮"),
    ("CCCCCCC(C)=O", "octan-2-one", None),
    ("CCCCCCCCC(C)=O", "decan-2-one", None),
    # negative: alkane / acid / alcohol must not become ketones
    ("CCC", "propane", "丙烷"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_ketone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
