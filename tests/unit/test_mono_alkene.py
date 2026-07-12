# IUPAC: P-31.1
# Layer: L1,L2,L4,L5
"""Simple acyclic monoalkenes: ane→ene; lowest double-bond locant.

C2/C3 omit locant (ethene/propene); C≥4 use stem-loc-ene / 位次-烯.
Parent chain must include the unique non-aromatic C=C.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic monoalkenes
    ("C=C", "ethene", "乙烯"),
    ("C=CC", "propene", "丙烯"),
    ("C=CCC", "but-1-ene", "丁-1-烯"),
    ("CC=CC", "but-2-ene", "丁-2-烯"),
    ("C=CCCC", "pent-1-ene", "戊-1-烯"),
    ("CC=CCC", "pent-2-ene", "戊-2-烯"),
    ("C=CCCCC", "hex-1-ene", "己-1-烯"),
    # negative: alkane / alcohol / aldehyde / ketone / acid must not become alkenes
    ("CCC", "propane", "丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_alkene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
