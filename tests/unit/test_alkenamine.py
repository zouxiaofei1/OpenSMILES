# IUPAC: P-62.2 / P-31.1
# Layer: L2,L4,L5
"""Unsaturated primary amines (alkenamines).

Amines use the segment-style ene insertion (same as alcohols):
but-3-en-1-amine. Previously the amine entry had no ene segment, so double
bonds were silently dropped (butan-1-amine).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # C≥3：烯/胺位次需区分位置异构 → 保留 (P-62.2.6.2)
    ("C=CCCN", "but-3-en-1-amine", "丁-3-烯-1-胺"),
    ("CC=CCN", "but-2-en-1-amine", "丁-2-烯-1-胺"),
    # C2：烯只能 1-2、胺后缀锚定 1，位次无歧义省略融合 (P-14.3.4)
    ("C=CN", "ethenamine", "乙烯胺"),
    # #502：短链省略仅去掉两个 '1'，另一端硝基位次 2 仍保留
    (r"C(C1=CC=CC=C1)N/C=C\[N+](=O)[O-]",
     "(1Z)-N-benzyl-2-nitroethenamine", "(1Z)-N-苄基-2-硝基乙烯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenamine(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
