# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4,L5
"""End-to-end benzothiazole / benzoxazole naming (retained parent P-22.2.1 / P-25).

Recognition now goes through the retained-template / ring_core path; the former
data-driven `_try_di13_fused56` producer engine was removed as dead code.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_E2E = [
    ("c1nc2ccccc2s1", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("c1ccc2ocnc2c1", "1,3-benzoxazole", "1,3-苯并噁唑"),
    ("Nc1nc2ccccc2s1", "1,3-benzothiazol-2-amine", "2-氨基苯并噻唑"),
    ("Brc1ccc2ncoc2c1", "6-bromo-1,3-benzoxazole", "6-溴-1,3-苯并噁唑"),
    ("c1coc2ccccc12", "benzofuran", "苯并呋喃"),
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_e2e_fused56_di13(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
