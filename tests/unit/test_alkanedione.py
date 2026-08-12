# IUPAC: P-64.2.1
# Layer: L2,L4,L5
"""Unsubstituted open-chain saturated alkanediones (exactly two ketones).

Parent chain through both ketone carbons; name alkane-{a},{b}-dione /
{首字}-{a},{b}-二酮 (mono-ketone Chinese style).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted open-chain diones
    ("CC(=O)CC(=O)C", "pentane-2,4-dione", "戊-2,4-二酮"),
    ("CC(=O)C(C)=O", "butane-2,3-dione", "丁-2,3-二酮"),
    ("CC(=O)CCC(=O)C", "hexane-2,5-dione", "己-2,5-二酮"),
    ("CC(=O)CCCC(=O)C", "heptane-2,6-dione", "庚-2,6-二酮"),
    ("CC(=O)CCCCC(=O)C", "octane-2,7-dione", "辛-2,7-二酮"),
    ("CCC(=O)CC(=O)CC", "heptane-3,5-dione", "庚-3,5-二酮"),
    # negative: monoketones / cycloketone must not break
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("CCCC(=O)C", "pentan-2-one", "戊-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanedione(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
