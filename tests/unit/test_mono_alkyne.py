# IUPAC: P-31.1
# Layer: L1,L2,L4,L5
"""Simple acyclic monoalkynes: ane→yne; lowest triple-bond locant.

C2 retained acetylene/乙炔; C3 propyne/丙炔 (omit locant);
C≥4 use stem-loc-yne / 位次-炔. Parent chain includes unique C≡C.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic monoalkynes
    ("C#CC", "propyne", "丙炔"),
    ("C#CCC", "but-1-yne", "丁-1-炔"),
    ("CC#CC", "but-2-yne", "丁-2-炔"),
    ("C#CCCC", "pent-1-yne", "戊-1-炔"),
    ("CC#CCC", "pent-2-yne", "戊-2-炔"),
    ("C#CCCCC", "hex-1-yne", "己-1-炔"),
    # negative: alkane / alkene / alcohol / aldehyde / ketone / acid
    ("CCC", "propane", "丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_alkyne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
