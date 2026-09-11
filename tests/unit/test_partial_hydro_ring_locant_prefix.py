# IUPAC: P-22.2.2, P-25.1
# Layer: L2
"""局部不饱和保留名的词干已含加氢前缀时，locant 前缀须留在组分名前，不得前移。

`_TEMPLATES` 里 dihydrothiazole 的词干是 `4,5-dihydro-1,3-thiazole`：`1,3-` 必须
紧贴 `thiazole`。若按「词干为裸名、locant 前缀在词首注入」的路径再补一次，就得到
`1,3-4,5-dihydro-1,3-thiazole`（位次重复）。同族 dihydropyrrole / dihydroimidazole
的 `1H-` 同样嵌在词中，靠 `1H-` 只在环含未取代芳香 NH 时才注入而侥幸未暴露。

对照条守住边界：词干不含加氢前缀者（thiazolidine / dioxolane / benzothiazole）
仍须在词首注入 locant 前缀。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")
CASES = [
    # 用户报例：4,5-二氢-1,3-噻唑-4-羧酸
    (
        "O=C(O)[C@H]1CSC(c2c[nH]c3ccccc23)=N1",
        "(4S)-2-(1H-indol-3-yl)-4,5-dihydro-1,3-thiazole-4-carboxylic acid",
        "(4S)-2-(1H-吲哚-3-基)-4,5-二氢-1,3-噻唑-4-羧酸",
    ),
    # 裸环：词干原样
    ("C1=NCCS1", "4,5-dihydro-1,3-thiazole", "4,5-二氢-1,3-噻唑"),
    # 同族（1H- 嵌在词中）不得前移
    ("C1C=CCN1", "2,5-dihydro-1H-pyrrole", "2,5-二氢-1H-吡咯"),
    ("C1=NCCN1", "4,5-dihydro-1H-imidazole", "4,5-二氢-1H-咪唑"),
    # 对照：词干不含加氢前缀，locant 前缀仍在词首注入
    ("C1NCCS1", "1,3-thiazolidine", "1,3-噻唑烷"),
    ("C1COCO1", "1,3-dioxolane", "1,3-二氧戊环"),
    ("c1ccc2scnc2c1", "1,3-benzothiazole", "1,3-苯并噻唑"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_partial_hydro_ring_locant_prefix_not_duplicated(smiles: str, en: str, zh: str) -> None:
    """部分不饱和保留名的 locant 前缀只出现一次，且留在组分名前。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
