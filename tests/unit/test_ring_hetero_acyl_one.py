# IUPAC: P-65.6.3.5
# Layer: L1
"""环内杂原子酰基按环酮命名：-C(=O)-O- / -C(=O)-S- 在环内（内酯/硫代内酯）同作 -one。

内酯（环内 O）依 P-65.6.3.5 按杂环 -one 命名，环内 O 换 S 同理；但 L1 的
`_is_lactone_carbon` 只认 O，硫代内酯的羰基不被识别为任何官能团，母体降解为
alkane、羰基氧无人归属（coverage 不完整），经 no_coverage_gate 兜底后输出丢失
羰基氧的稠合母体氢化物名。本轮把「羰基碳的环内杂原子邻居（N/O/S）」收敛为
`_has_ring_hetero_neighbor`，使环内 S 与环内 O/N 同样作环酮。

末三条为 n_c==0 兄弟分支（环碳酸酯/环氨基甲酸酯）的守卫：该分支同样依赖内酯
判定，改动 `_is_ketone_carbon` 时若漏掉其 `lactone` 依赖会整体抛 NameError。

pos2/pos3 的 S 稠合苯并母体现取通用稠合名 benzo[b]thian-*（ring_scaffold 尚未登记
thiochromene/thiochromane 保留母体），与 O 版 chromen-2-one 的词干选择不同，此处按
当前实现断言。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: 5 元环内硫代内酯（isobenzothiophenone），与 3H-2-benzofuran-1-one 同构
    ("C1C2C=CC=CC=2C(=O)S1", "3H-benzo[c]thiophen-1-one", "3H-苯并[c]噻吩-1-酮"),
    # positive: 6 元不饱和硫代内酯（苯并稠合）与饱和硫代内酯
    ("O=C1C=Cc2ccccc2S1", "benzo[b]thian-2-one", "苯并[b]四氢噻喃-2-酮"),
    ("O=C1CCc2ccccc2S1", "3,4-dihydrobenzo[b]thian-2-one", "3,4-二氢苯并[b]四氢噻喃-2-酮"),
    # positive: 单环饱和 6 元硫代内酯，与 O 版 oxan-2-one 逐字对应
    ("O=C1SCCCC1", "thian-2-one", "四氢噻喃-2-酮"),
    # negative: 环内 S 但无羰基 —— 不得改判为酮
    ("C1CCSCC1", "thiane", "四氢噻喃"),
    # negative: 非环硫醚 —— 环内 S 分支不得误伤
    ("CSC", "methylsulfanylmethane", "甲硫基甲烷"),
    # negative: O 版内酯与 N 版内酰胺的既有命名不得被本轮改动带偏
    ("O=C1OCCCC1", "oxan-2-one", "氧杂环己烷-2-酮"),
    ("C1C2C=CC=CC=2C(=O)N1", "2,3-dihydroisoindol-1-one", "2,3-二氢异吲哚-1-酮"),
    # negative guard: n_c==0 环碳酸酯 / 环氨基甲酸酯 / 环脲保持原命名（兄弟分支不得被带崩）
    ("O=C1OCCO1", "1,3-dioxolan-2-one", "1,3-二氧戊环-2-酮"),
    ("O=C1OCCN1", "1,3-oxazolidin-2-one", "1,3-噁唑烷-2-酮"),
    ("O=C1NCCN1", "imidazolidin-2-one", "咪唑烷-2-酮"),
    (
        "Cl[C@@H]1OC(O[C@H]1C(Cl)(Cl)Cl)=O",
        "(4S,5R)-4-chloro-5-(trichloromethyl)-1,3-dioxolan-2-one",
        "(4S,5R)-4-氯-5-(三氯甲基)-1,3-二氧戊环-2-酮",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ring_hetero_acyl_one(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
