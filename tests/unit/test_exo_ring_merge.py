# IUPAC: P-65.2.2 / P-66.6.1.1.3 / P-65.1.7.2
# Layer: L5
#
# 环外主基（-COOH/-CHO/-COOR/-CONH2/-CN/-CO-）统一由链引擎产出：环母体词干 + -carb- 后缀。
#
# 环外后缀拼的是**完整母体氢化物**（cyclohexane + carboxylic acid），环内不饱和落进词干并按段式给位次
# （cyclohex-2-yne-1-carboxylic acid / 环己-2-炔-1-羧酸），与链引擎的 ene/yne 段共用一条管线——
# 旧 exocyclic worker 只读 ene_locants 且自带词干拼装，炔键与位次省略都漏。
# 附：环上另有被编号前缀时主基位次不可省略（2-methylcyclohexane-1-carboxylic acid）。
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# 新增：环内炔键此前被整段丢弃（旧 worker 只读 ene_locants）
EXO_RING_UNSAT = [
    ("C1C(C(=O)O)C#CCC1", "cyclohex-2-yne-1-carboxylic acid", "环己-2-炔-1-羧酸"),
]

# 回归：环烯/环烷/稠环/杂环/多取代环外主基的既有正确形态
EXO_RING_KEEP = [
    ("C1CCC(C(=O)O)CC1", "cyclohexanecarboxylic acid", "环己烷羧酸"),
    ("O=C(O)C1CC1", "cyclopropanecarboxylic acid", "环丙烷羧酸"),
    ("C1CC(C(=O)O)CC=C1", "cyclohex-3-ene-1-carboxylic acid", "环己-3-烯-1-羧酸"),
    # 单烯 1 位：EN 省位次、中文保留（P-31.1.2；金标 cyclohexene-1-carboxylic acid / 环己-1-烯-1-羧酸）
    ("C1CC=C(C(=O)O)CC1", "cyclohexene-1-carboxylic acid", "环己-1-烯-1-羧酸"),
    ("C=C(C)C1CC=C(C=O)CC1", "4-prop-1-en-2-ylcyclohexene-1-carbaldehyde",
     "4-丙-1-烯-2-基环己-1-烯-1-甲醛"),
    ("O=C([O-])C1=CC(=O)[C@@H](O)[C@H](O)C1",
     "(4S,5R)-4,5-dihydroxy-3-oxocyclohexene-1-carboxylate",
     "(4S,5R)-4,5-二羟基-3-氧代环己-1-烯-1-羧酸根"),
    ("CC1CCCCC1C(=O)O", "2-methylcyclohexane-1-carboxylic acid", "2-甲基环己烷-1-羧酸"),
    ("O=CC1CCCCC1", "cyclohexanecarbaldehyde", "环己烷甲醛"),
    ("O=CC1CC=CCC1", "cyclohex-3-ene-1-carbaldehyde", "环己-3-烯-1-甲醛"),
    ("NC(=O)C1CCCCC1", "cyclohexanecarboxamide", "环己烷甲酰胺"),
    ("CNC(=O)C1CCCCC1", "N-methylcyclohexanecarboxamide", "N-甲基环己烷甲酰胺"),
    ("N#CC1CCCCC1", "cyclohexanecarbonitrile", "环己烷甲腈"),
    ("COC(=O)C1CCCCC1", "methyl cyclohexanecarboxylate", "环己烷羧酸甲酯"),
    ("O=C(O)C1CCC(C(=O)O)CC1", "cyclohexane-1,4-dicarboxylic acid", "环己烷-1,4-二羧酸"),
    # 稠环（scaffold_id 为 carbocycle 但带 fused_tree）：母体取保留稠环名，不得退化成环烷烃（中文稠环词干另有既存偏差，只守 EN）
    ("Cc1c(C(=O)O)c(O)cc2cc3c(c(O)c12)C(=O)c1c(O)cc(O)cc1C3=O",
     "3,8,10,12-tetrahydroxy-1-methyl-6,11-dioxotetracene-2-carboxylic acid", None),
    ("N1C=C(C2=CC=C(C=C12)C=O)C=O", "1H-indole-3,6-dicarbaldehyde", "1H-吲哚-3,6-二甲醛"),
    ("O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1",
     "benzene-1,3,5-tricarboxylic acid", "苯-1,3,5-三羧酸"),
]

# 近邻负例：保留名（苯/杂环/稠环单取代）、开链酸 —— 不得被环外通用式抢占
EXO_RING_NEIGHBOR = [
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("O=C(O)c1ccco1", "furan-2-carboxylic acid", "呋喃-2-羧酸"),
    ("O=C(O)c1ccc2ccccc2c1", "naphthalene-2-carboxylic acid", "萘-2-羧酸"),
    ("OC(=O)CCCCCC", "heptanoic acid", "庚酸"),
]


@pytest.mark.parametrize("smiles,en,zh", EXO_RING_UNSAT)
def test_exo_ring_unsat_locants(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


@pytest.mark.parametrize("smiles,en,zh", EXO_RING_KEEP)
def test_exo_ring_keep(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


@pytest.mark.parametrize("smiles,en,zh", EXO_RING_NEIGHBOR)
def test_exo_ring_neighbor(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh
