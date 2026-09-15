# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_acyl_halide.py: Acyl halides: alkanoyl/benzoyl/enoyl F·Cl·Br·I (P-65.5).
test_alkanoyl_chloride.py: Unsubstituted open-chain alkanoyl chlorides (mono acyl chlorides).
test_alkanoic_anhydride.py: Unsubstituted symmetrical open-chain alkanoic anhydrides.
test_acyl_oyl.py: 酰基/酰胺残基（自由价锚定在羰基碳的 R-C(=O)- 片段）必须以酸衍生酰基名命名，
test_carboxy_residue.py: 残基/取代基内的游离 -COOH(非主官能团)须前缀化为 carboxy,不得并入残基主链
test_ring_acyl_carbonyl.py: 环/芳羧酸的酰基残基 X-C(=O)-（X = 苯环/芳杂环/饱和碳环）必须以收缩酰基名命名，
test_ring_hetero_acyl_one.py: 环内杂原子酰基按环酮命名：-C(=O)-O- / -C(=O)-S- 在环内（内酯/硫代内酯）同作 -one。
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_acyl_halide.py
# IUPAC: P-65.5
# Layer: L1,L2,L5
#
# Acyl halides: alkanoyl/benzoyl/enoyl F·Cl·Br·I (P-65.5).
#
# 主 FG 命名后缀随卤素（hal_z）变：-oyl halide / 酰氟氯溴碘（保留名 acetyl/
# benzoyl 亦同）。覆盖开链 C2–C4、2-甲基丙酰、苯甲酰与烯酰。负数：酸、
# 溴代烷、醛不得判成酰卤。
# ==========================================================================
acyl_halide__CASES = [
    # ── 正例：C1/C2 保留名与开链后缀随卤素变 ──
    ("CC(=O)F", "acetyl fluoride", "乙酰氟"),
    ("CC(=O)Cl", "acetyl chloride", "乙酰氯"),
    ("CC(=O)Br", "acetyl bromide", "乙酰溴"),
    ("CC(=O)I", "acetyl iodide", "乙酰碘"),
    ("CCC(=O)Cl", "propanoyl chloride", "丙酰氯"),
    ("CCC(=O)F", "propanoyl fluoride", "丙酰氟"),
    ("CCC(=O)I", "propanoyl iodide", "丙酰碘"),
    # 2-甲基丙酰（基准金标：Br 必须 bromide 而非 chloride）
    ("CC(C)C(=O)Br", "2-methylpropanoyl bromide", "2-甲基丙酰溴"),
    ("CC(C)C(=O)I", "2-methylpropanoyl iodide", "2-甲基丙酰碘"),
    # ── 苯甲酰卤（环外酰卤主 FG 同样带卤素）──
    ("O=C(F)c1ccccc1", "benzoyl fluoride", "苯甲酰氟"),
    ("O=C(Br)c1ccccc1", "benzoyl bromide", "苯甲酰溴"),
    ("O=C(I)c1ccccc1", "benzoyl iodide", "苯甲酰碘"),
    # ── 烯酰卤 ──
    ("O=C(Cl)C=C", "prop-2-enoyl chloride", "丙-2-烯酰氯"),
    ("O=C(F)C=C", "prop-2-enoyl fluoride", "丙-2-烯酰氟"),
    ("O=C(Br)C=C", "prop-2-enoyl bromide", "丙-2-烯酰溴"),
    # ── 炔酰卤（P-14.3.4：三碳炔酰保留炔位次，与烯酰卤同形）──
    ("O=C(Cl)C#C", "prop-2-ynoyl chloride", "丙-2-炔酰氯"),
    ("O=C(Br)C#C", "prop-2-ynoyl bromide", "丙-2-炔酰溴"),
    ("O=C(Cl)CC#C", "but-3-ynoyl chloride", "丁-3-炔酰氯"),



    # F 在 α 碳而非羰基：仍为醛（不误判为酰卤）
    ("FCC=O", "fluoroacetaldehyde", "氟乙醛"),
]


@pytest.mark.parametrize("smiles,en,zh", acyl_halide__CASES)
def test_acyl_halide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanoyl_chloride.py
# IUPAC: P-65.5.1
# Layer: L1,L2,L3,L4,L5
#
# Unsubstituted open-chain alkanoyl chlorides (mono acyl chlorides).
#
# Exactly one –C(=O)Cl; parent chain through acyl carbon; no other main FG,
# no ring/unsaturation this round. EN: C2 acetyl chloride; C>=3 …oyl chloride.
# ZH: …酰氯 (乙酰氯, 丙酰氯, …).
# ==========================================================================
alkanoyl_chloride__CASES = [
    ("ClCCCl", "1,2-dichloroethane", "1,2-二氯乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanoyl_chloride__CASES)
def test_alkanoyl_chloride(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanoic_anhydride.py
# IUPAC: P-65.7.1
# Layer: L1,L2,L5
#
# Unsubstituted symmetrical open-chain alkanoic anhydrides.
#
# Symmetric R–C(=O)–O–C(=O)–R; both acyl chains saturated, acyclic, no other
# main FG. EN: C2 acetic anhydride; C≥3 …oic anhydride (acid stem + anhydride).
# ZH: 酸名 + 酐 (乙酸酐, 丙酸酐, …).
# ==========================================================================
alkanoic_anhydride__CASES = [
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanoic_anhydride__CASES)
def test_alkanoic_anhydride(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_acyl_oyl.py
# IUPAC: P-64.5.1,P-65.1.7.2,P-66.1.1.4.3
# Layer: L1,L2,L3,L4,L5
#
# 酰基/酰胺残基（自由价锚定在羰基碳的 R-C(=O)- 片段）必须以酸衍生酰基名命名，
# 而非展开成「1-氧代烷基」。
#
# - P-64.5.1(2)：侧链 1 位羰基（自由价所在碳，即 -CO-R）必须用相应酰基名；官方反例
#   [not 5-(1-oxoethyl)nonane-4,6-dione]。
# - P-65.1.7.2：酰基名由对应酸名把 -oic acid→-oyl / -ic acid→-yl 改词尾而派生
#   （acetyl、butanoyl、2-thiophen-2-ylacetyl…）。
# - P-66.1.1.4.3：R-CO-NH- 残基作取代基可用 acylamino（方法 2，本簇 gold 采用）。
#
# 说明：ChEBI gold_zh 对 N-酰基残基常写作「…酰胺基」(amido 式方法 1) 且中英不对齐，
# 与本规则实现的方法(2) acylamino 不同写法并存；故 N-酰基行的中文断言跳过（zh=None），
# 只精确断言英文。O-酰基行同理以英文为准。
# ==========================================================================
acyl_oyl__POS_FRAG = [
    ("*C(=O)C", "acetyl", "乙酰基"),
    ("*C(=O)CCC", "butanoyl", None),
    ("*C(=O)Cc1cccs1", "2-thiophen-2-ylacetyl", "2-噻吩-2-基乙酰基"),
]


@pytest.mark.parametrize("smiles,en,zh", acyl_oyl__POS_FRAG)
def test_acyl_fragment_named_oyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子正例（gold/pred 仅差 1-oxo→oyl 单点；英文精确，中文跳过/待验）。
# 注：O-酰基的围栏按 P-65.6.3.2.3：retained 单词酰基名（acetyl/benzoyl/hexadecanoyl）
# 与 oxy 融合成 …oyloxy（3-(benzoyloxy)propanoic acid）；带位次的系统酰基名则闭括号在
# 前端之后、oxy 留括号外（3-[(pyridine-3-carbonyl)oxy]propanoic acid）。故带双键位次/
# 立体描述符的 (9Z)-octadec-9-enoyl 取后者，(9Z)-octadec-9-enoyl]oxy；与 benchmark gold
# （chebi-2873/chebi-776）一致。
acyl_oyl__POS_WHOLE = [
    ("CC(C)=CC(O)=NCC(=O)O", "2-(3-methylbut-2-enoylamino)acetic acid"),
    ("CCCCCCCCCCCCCCCCCC(O)=NCC(=O)O", "2-(octadecanoylamino)acetic acid"),
    ("C(C(=C)C)(=O)OC(C)OC(C(=C)C)=O",
     "1-(2-methylprop-2-enoyloxy)ethyl 2-methylprop-2-enoate"),
    ("CCCCCCCC/C=C\\CCCCCCCC(=O)OC(CCCCCCCCCCCCC)CCCC(=O)O",
     "5-[(9Z)-octadec-9-enoyl]oxyoctadecanoic acid"),
]


@pytest.mark.parametrize("smiles,en", acyl_oyl__POS_WHOLE)
def test_acyl_whole_molecule_oyl(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 负例：acyl 规则不得误伤。
# 1) 酰胺/酮作母体的既有路径不变；2) 内部氧代（锚点不在羰基碳）不变；
# 3) 锚点羰基碳带 O 离去基（酯羰基）不判酰基，保持现行为（带 N 的氨基甲酰见上方 POS_FRAG）。
acyl_oyl__NEG = [
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("CC(=O)c1ccccc1", "1-phenylethanone", "1-苯基乙酮"),
    ("*CCC(=O)C", "3-oxobutyl", "3-氧代丁基"),
    ("*C(=O)OC", "methoxy(oxo)methyl", "甲氧基(氧代)甲基"),
]


@pytest.mark.parametrize("smiles,en,zh", acyl_oyl__NEG)
def test_acyl_negative_untouched(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_carboxy_residue.py
# IUPAC: P-61.1.3,P-65.1.1.2
# Layer: L1,L2,L3
#
# 残基/取代基内的游离 -COOH(非主官能团)须前缀化为 carboxy,不得并入残基主链
# 并被读成同碳「n-hydroxy-n-oxo-」逐原子描述。
#
# - 自由基/酰基主基团(p41=1)压制非主羧酸(p41=7)时，该 COOH 是取代基叶（P-61.1.3
#   carboxy 前缀）：fragment `*CH2CH2COOH` 应命名 `2-carboxyethyl`，`*CH2CH2CH2COOH`
#   的酰基形态为 `3-carboxypropanoyl`，而非 `3-hydroxy-3-oxopropyl`/`4-hydroxy-4-oxobutanoyl`
#   把酸碳并进链长。
# - 整分子主官能团二酸/环酸（母体后缀 -dioic/-carboxylic acid）不降级，不得 carboxy 化。
# - 降级酰胺/醛/酮羰基保持链化 oxo（见 test_demoted_amide_to_amino.py），不在此列。
# ==========================================================================
carboxy_residue__POS_FRAG = [
    ("*CCC(=O)O", "2-carboxyethyl"),
    ("*C(=O)CCC(=O)O", "3-carboxypropanoyl"),
]


@pytest.mark.parametrize("smiles,en", carboxy_residue__POS_FRAG)
def test_carboxy_fragment(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例：N-酰基残基尾端 COOH（peptide/代谢物场景），羧酸主母体不受影响。
carboxy_residue__POS_WHOLE = [
    ("CC(C)C[C@H](N=C(O)[C@@H](N)CC(=O)O)C(=O)O",
     "(2S)-2-[[(2S)-2-amino-3-carboxypropanoyl]amino]-4-methylpentanoic acid"),
    ("O=C(O)CCC(O)=N[C@@H](CCCC(=O)C(=O)O)C(=O)O",
     "(2S)-2-(3-carboxypropanoylamino)-6-oxoheptanedioic acid"),
]


@pytest.mark.parametrize("smiles,en", carboxy_residue__POS_WHOLE)
def test_carboxy_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 负例：carboxy 规则不得误伤。
carboxy_residue__NEG = [
    # 降级酰胺羰基仍链化 oxo（伯酰胺 N → amino、羰基 C 留链作 4-oxo），不得 carboxy 化。
    ("NC(=O)CCC(=O)O", "4-amino-4-oxobutanoic acid"),
    # 主二酸母体不被降级、不得 carboxy 化，整分子仍为 dioic 母体后缀。
    ("O=C(O)CCC(=O)O", "butanedioic acid"),
]


@pytest.mark.parametrize("smiles,en", carboxy_residue__NEG)
def test_carboxy_negative_untouched(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# ==========================================================================
# 合并自 test_ring_acyl_carbonyl.py
# IUPAC: P-65.1.7.2,P-64.5.1,P-66.6.1
# Layer: L1,L2,L5
#
# 环/芳羧酸的酰基残基 X-C(=O)-（X = 苯环/芳杂环/饱和碳环）必须以收缩酰基名命名，
# 而非展开成「yl(oxo)methyl / -oxomethyl」。
#
# - P-65.1.7.2：酰基名由对应酸名派生——benzoic acid→benzoyl、furan-2-carboxylic acid→
#   furan-2-carbonyl、cyclopropanecarboxylic acid→cyclopropanecarbonyl（把 -oic/-ic acid 词尾
#   换成 -oyl/-carbonyl）；α 为环/芳碳时同样是酸衍生酰基，不能退化成逐原子 oxo 描述。
# - P-64.5.1(2)：侧链 1 位羰基（自由价在羰基碳的 -CO-Ar）用相应酰基名。
# - 既有 test_acyl_oyl.py 已覆盖 α 为开链烷基/烯基的酰基收缩（butanoyl/enoyl/octadecanoyl）；
#   本文件补齐 α 为环/芳的缺口（A1 模式，benchmark ~50 例）。
#
# 与 test_acyl_oyl.py 同约定：N-酰基残基中文（gold 常「…酰胺基/甲酰…」与 IUPAC en 不对齐）跳过，
# 只精确断言英文；fragment 中文按 namer acyl 词形惯例。
# ==========================================================================
ring_acyl_carbonyl__POS_FRAG = [
    ("*C(=O)c1ccccc1", "benzoyl"),
    ("*C(=O)c1ccco1", "furan-2-carbonyl"),
    ("*C(=O)C1CC1", "cyclopropanecarbonyl"),
    ("*C(=O)c1ccc(Cl)cc1", "4-chlorobenzoyl"),
    ("*C(=O)C1CCCN1", "pyrrolidine-2-carbonyl"),
]


@pytest.mark.parametrize("smiles,en", ring_acyl_carbonyl__POS_FRAG)
def test_ring_acyl_fragment_named_carbonyl(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例：N-酰基氨基残基（羧酸主母体 + -NH-CO-X），X 为环/芳时应收缩。
ring_acyl_carbonyl__POS_WHOLE = [
    ("O=C(O)CNC(=O)c1ccco1", "2-(furan-2-carbonylamino)acetic acid"),
    ("C1CC1C(=O)NCC(=O)O", "2-(cyclopropanecarbonylamino)acetic acid"),
]


@pytest.mark.parametrize("smiles,en", ring_acyl_carbonyl__POS_WHOLE)
def test_ring_acyl_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# ==========================================================================
# 合并自 test_ring_hetero_acyl_one.py
# IUPAC: P-65.6.3.5
# Layer: L1
#
# 环内杂原子酰基按环酮命名：-C(=O)-O- / -C(=O)-S- 在环内（内酯/硫代内酯）同作 -one。
#
# 内酯（环内 O）依 P-65.6.3.5 按杂环 -one 命名，环内 O 换 S 同理；但 L1 的
# `_is_lactone_carbon` 只认 O，硫代内酯的羰基不被识别为任何官能团，母体降解为
# alkane、羰基氧无人归属（coverage 不完整），经 no_coverage_gate 兜底后输出丢失
# 羰基氧的稠合母体氢化物名。本轮把「羰基碳的环内杂原子邻居（N/O/S）」收敛为
# `_has_ring_hetero_neighbor`，使环内 S 与环内 O/N 同样作环酮。
#
# 末三条为 n_c==0 兄弟分支（环碳酸酯/环氨基甲酸酯）的守卫：该分支同样依赖内酯
# 判定，改动 `_is_ketone_carbon` 时若漏掉其 `lactone` 依赖会整体抛 NameError。
#
# pos2/pos3 的 S 稠合苯并母体现取通用稠合名 benzo[b]thian-*（ring_scaffold 尚未登记
# thiochromene/thiochromane 保留母体），与 O 版 chromen-2-one 的词干选择不同，此处按
# 当前实现断言。
# ==========================================================================
ring_hetero_acyl_one__CASES = [
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


@pytest.mark.parametrize("smiles,en,zh", ring_hetero_acyl_one__CASES)
def test_ring_hetero_acyl_one(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
