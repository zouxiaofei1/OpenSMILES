# IUPAC: P-62.2 / P-15.4.1
# Layer: L5
"""桥后缀（氨基/氧基/硫基）前中文烃基名的「基」字去留。

规则：桥后缀直接相连的烃基名省略尾「基」字（甲基→甲、环己基→环己、
丙-2-基→丙-2-），与 tiers gold 口径一致。母体为胺时 N-取代前缀
（N,N-二甲基苯胺）不受此规则影响，仍保留「基」。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")
CASES = [
    # ── 正例：N-多取代简单烷氨基作为取代基，中文省「基」 ──
    ("CN(C)CCO", "2-(dimethylamino)ethanol", "2-(二甲氨基)乙醇"),
    ("CCN(CC)c1ccccc1C(=O)O", "2-(diethylamino)benzoic acid", "2-(二乙氨基)苯甲酸"),
    ("CCCN(CCC)c1ccccc1C(=O)O", "2-(dipropylamino)benzoic acid", "2-(二丙氨基)苯甲酸"),
    # ── 近邻：位次化/环系烃基名同样省「基」（tiers gold：丙-2-氧基、环己氧基、环丙硫基） ──
    ("CC(C)NCCO", "2-(propan-2-ylamino)ethanol", "2-(丙-2-氨基)乙醇"),
    ("C1CCCCC1NCCO", "2-(cyclohexylamino)ethanol", "2-(环己氨基)乙醇"),
    # ── 近邻负例：胺为母体时的 N-取代前缀仍保留「基」 ──
    ("CN(C)c1ccccc1", "N,N-dimethylaniline", "N,N-二甲基苯胺"),
    ("CN(C)CCN(C)C", "N,N,N',N'-tetramethylethane-1,2-diamine",
     "N,N,N',N'-四甲基乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_zh_bridge_root(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
