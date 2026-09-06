# IUPAC: P-31.1.1.1 / P-31.1.1.2
# Layer: L2,L4,L5
"""混合烯炔与多炔命名：ene 前 yne 后、词干插 a、多重键集合编号优先。

纯烃（alkane）、FG 段式（醇）与融合式（酸）三路均覆盖；多烯+多炔、带侧链一并。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # 混合烯炔（纯烃）：ene 前 yne 后，双键得低位
    ("C#CCC=C", "pent-1-en-4-yne", "戊-1-烯-4-炔"),
    ("C=CC#C", "but-1-en-3-yne", "丁-1-烯-3-炔"),
    # 3 位侧基 =CH2 是外环亚甲基（P-56.4 methylidene）；旧 'methyl' 丢双键致式量错
    # （C=C(C#C)C=C 为 C6H6，饱和甲基对应 C6H10，不同分子）。
    ("C=C(C#C)C=C", "3-methylidenepent-1-en-4-yne", "3-亚甲基戊-1-烯-4-炔"),
    # 多炔：词干插 a，MULT 后缀
    ("C#CCCC#C", "hexa-1,5-diyne", "己-1,5-二炔"),
    ("C#CC#C", "buta-1,3-diyne", "丁-1,3-二炔"),
    ("C#CCC#CC", "hexa-1,4-diyne", "己-1,4-二炔"),
    # 多烯+多炔混合：dien + diyne 组合段
    ("C=C=CC#CC#C", "hepta-1,2-dien-4,6-diyne", "庚-1,2-二烯-4,6-二炔"),
    # FG 段式（醇）
    ("C#CCC(O)C=C", "hex-1-en-5-yn-3-ol", "己-1-烯-5-炔-3-醇"),
    # FG 融合式（酸）
    ("C=CC#CCC(=O)O", "hex-5-en-3-ynoic acid", "己-5-烯-3-炔酸"),
    # 多 FG（多醇/多胺）：词干走主路径，多炔/混合一并支持
    ("OCC=CCO", "but-2-ene-1,4-diol", "丁烷-2-烯-1,4-二醇"),
    ("OCC#CC#CCO", "hexa-2,4-diyne-1,6-diol", "己烷-2,4-二炔-1,6-二醇"),
    ("NCC=CCN", "but-2-ene-1,4-diamine", "丁烷-2-烯-1,4-二胺"),
    ("NCCC#CCN", "pent-2-yne-1,5-diamine", "戊烷-2-炔-1,5-二胺"),
    ("NC=CC#CCN", "pent-1-en-3-yne-1,5-diamine", "戊烷-1-烯-3-炔-1,5-二胺"),
    ("NCC=CC#CCN", "hex-2-en-4-yne-1,6-diamine", "己烷-2-烯-4-炔-1,6-二胺"),
    # 多 FG（多硫醇）：P-57 保留 e
    ("SCC=CCS", "but-2-ene-1,4-dithiol", "丁烷-2-烯-1,4-二硫醇"),
    ("SCCC#CCS", "pent-2-yne-1,5-dithiol", "戊烷-2-炔-1,5-二硫醇"),
    ("SC=CC#CCS", "pent-1-en-3-yne-1,5-dithiol", "戊烷-1-烯-3-炔-1,5-二硫醇"),
    # 对照 negative：纯烯/纯炔不得走混合路径
    ("C#C", "ethyne", "乙炔"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_enyne_polyyne(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
