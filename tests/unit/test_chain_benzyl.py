# IUPAC: P-29.3 / P-63.1 / P-65.1 / P-62.2
# Layer: L2,L3,L5
"""Benzyl / phenyl on chain acid / amine parents (plan 3.4)."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # diacid + C-benzyl / phenyl (open chain no longer walks into arene)
    (
        "O=C(O)CC(Cc1ccccc1)C(=O)O",
        "2-benzylbutanedioic acid",
        "2-苄基丁二酸",
    ),
    # phenylmethanamine system name；乙胺 C1 有可取代 H，2- 不可省略（P-14.3.4.4）
    ("c1ccc(CCN)cc1", "2-phenylethanamine", "2-苯基乙胺"),
    # regressions
    ("CC(=O)O", "acetic acid", None),
    # P-14.3.4.4 异构体测试：乙醇 C1 有可取代 H，2- 不可省略
    ("NCCO", "2-aminoethanol", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_chain_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
