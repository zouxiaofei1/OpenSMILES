# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4,L5
"""End-to-end benzofuran / benzothiophene naming (retained parent P-22.2.1 / P-25).

Recognition now goes through the retained-template / ring_core path; the former
data-driven `_try_mono_fused56` producer engine was removed as dead code.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

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
