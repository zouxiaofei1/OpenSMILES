# 合并自 9 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mult_prefix.py: Multiplicative prefixes penta–deca for identical substituents (P-14.2.1 / P-16.3).
test_mult_prefix_bis_di.py: 相同复合取代基重复时倍增前缀须用 bis/tris 而非 di/tri。
test_nprefix_locant_mixed.py: N 型取代基与 C 型取代基同词干混合组：N 成员 locant 必须渲染为 N（而非 0/羰基碳位/环连接位假数字）。
test_bridge_suffix_parens.py: O/S/N 桥取代基前缀的括注边界与嵌套层级。
test_simple_alkyl_prefix_paren.py: Linear n-alkyl as substituent prefix: NO parentheses when simple & locant-free.
test_omit_locant_c2_fg.py: C2 单取代母体 locant 省略收窄：醇/硫醇/胺 FG 端碳有可取代 H，2- 不可省略。
test_tert_pentyl.py: 2-methylbutan-2-yl (tert-pentyl / 1,1-dimethylpropyl) branched alkyl.
test_branched_alkyl_mode.py: Verify general vs pin mode name switching for branched alkyls.
test_monoalkyl_chain.py: Straight-chain monoalkyl (C1–C4) side-chain prefixes on alkane/alcohol parents.
"""
from __future__ import annotations

import pytest

from namepredict.layer5.assembler_prefixes import _front_needs_enclosure, _split_bridge_suffix
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import alkyl_alpha_key, normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mult_prefix.py
# IUPAC: P-14.2.1
# Layer: L5
#
# Multiplicative prefixes penta–deca for identical substituents (P-14.2.1 / P-16.3).
#
# Affiliated: P-61.3.1 halogen stems; existing di/tri/tetra must not regress.
# ==========================================================================
mult_prefix__CASES = [
    # positive: count ≥5 needs penta/hexa/… multiplicative prefixes
    # zh=None when gold has no Chinese (merged_benchmark empty chinese_name)
    ("ClC(Cl)C(Cl)(Cl)Cl", "1,1,1,2,2-pentachloroethane", None),
    ("ClC(Cl)(Cl)C(Cl)(Cl)Cl", "1,1,1,2,2,2-hexachloroethane", None),
    ("FC(F)(F)C(F)(F)F", "1,1,1,2,2,2-hexafluoroethane", None),
    ("FC(F)(F)C(F)(F)C(F)F", "1,1,1,2,2,3,3-heptafluoropropane", None),
    ("FC(F)(F)C(F)(F)C(F)(F)F", "1,1,1,2,2,3,3,3-octafluoropropane", None),
    # negative: di/tetra already correct — must not regress
    ("ClC(Cl)(Cl)Cl", "tetrachloromethane", "四氯甲烷"),
    ("BrCBr", "dibromomethane", None),
]


@pytest.mark.parametrize("smiles,en,zh", mult_prefix__CASES)
def test_mult_prefix(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_mult_prefix_bis_di.py
# IUPAC: P-16.3.2（乘法前缀选择：被取代/复合组分用 bis/tris，简单组分用 di/tri）
# Layer: L5 (assembler_prefixes)
#
# 相同复合取代基重复时倍增前缀须用 bis/tris 而非 di/tri。
#
# 修复前 namer 仅对含 'carboxy' 的词干用 bis，漏掉 hydroxymethyl（retained 组合叶）
# 与 recursive 复合前缀（trimethoxyphenyl/bromobutyl/hydroxyphenyl/hydroxyethoxy 等），
# 输出 di(...) 而 gold/IUPAC 为 bis(...)。以下 7 例修复后应与 gold 全串一致。
# ==========================================================================
mult_prefix_bis_di__FULL_CASES = [
    # A12 代表：retained 组合叶 hydroxymethyl 与 recursive 复合前缀
    ("C/C(=C\\CC/C(C)=C/CC/C=C(\\C)CC/C=C(\\C)CCC=C(CO)CO)CCC=C(CO)CO",
     "(6E,10E,14E,18E)-2,23-bis(hydroxymethyl)-6,10,15,19-tetramethyltetracosa-2,6,10,14,18,22-hexaene-1,24-diol"),
    ("OCC[C@@H]1C(CO)=C(CO)C[C@H]1O",
     "(1R,2R)-2-(2-hydroxyethyl)-3,4-bis(hydroxymethyl)cyclopent-3-en-1-ol"),
    ("COC1=C(C=CC(=C1OC)OC)\\C=C\\C(CC(\\C=C\\C1=C(C(=C(C=C1)OC)OC)OC)=O)=O",
     "(1E,6E)-1,7-bis(2,3,4-trimethoxyphenyl)hepta-1,6-diene-3,5-dione"),
    ("BrCCCCC1=CC=C(C=C1)CCCCBr",
     "1,4-bis(4-bromobutyl)benzene"),
    ("O=C(/C=C/c1ccc(O)cc1)CC(=O)/C=C/c1ccc(O)cc1",
     "(1E,6E)-1,7-bis(4-hydroxyphenyl)hepta-1,6-diene-3,5-dione"),
    ("CCCCCCCCCCCCCCCCCC(=O)OCCOCC(OCCO)[C@H]1OCC(OCCO)[C@H]1OCCO",
     "2-[2-[(2R,3R)-3,4-bis(2-hydroxyethoxy)oxolan-2-yl]-2-(2-hydroxyethoxy)ethoxy]ethyl octadecanoate"),
    # A12 家族第 7 条：方括号复杂基团 1-(2-methylbutan-2-yl)indol-3-yl 双取代
    ("CCC(C)(C)n1cc(C2=C(O)C(=O)C(c3cn(C(C)(C)CC)c4ccccc34)=C(O)C2=O)c2ccccc21",
     "2,5-dihydroxy-3,6-bis[1-(2-methylbutan-2-yl)indol-3-yl]cyclohexa-2,5-diene-1,4-dione"),
]


@pytest.mark.parametrize("smiles,en", mult_prefix_bis_di__FULL_CASES)
def test_bis_compound_mult(smiles: str, en: str) -> None:
    """复合取代基倍增应输出 bis(...) 且整串与 gold 一致。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# ==========================================================================
# 合并自 test_nprefix_locant_mixed.py
# IUPAC: P-45 / P-62 / P-66 (N- 与 C 数字位次并入同一乘数前缀)
# Layer: L5 (assembler_prefixes)
#
# N 型取代基与 C 型取代基同词干混合组：N 成员 locant 必须渲染为 N（而非 0/羰基碳位/环连接位假数字）。
#
# 触发：母体单胺/单酰胺，某词干同时出现在 C 位与该 N 上；合并组落入数字前缀通道，
# N 成员的假 locant（胺 N=0、酰胺 N 被 remap 到羰基碳/环连接位）被当普通数字打印。
# 期望形如 N,6,6-trimethyl / N,N,2-trimethyl / N,2-dimethyl / N,3-dimethyl（N 在最前）。
# ==========================================================================
nprefix_locant_mixed__FULL_CASES = [
    # 酰胺：C2-methyl + N,N-二甲基 → N,N,2-trimethyl（曾输出 1,1,2-）
    ("CC(C)C(=O)N(C)C", "N,N,2-trimethylpropanamide"),
    # 伯胺：C2-methyl + N,N-二甲基 → N,N,2-trimethyl（曾输出 0,0,2-）
    ("CC(C)CN(C)C", "N,N,2-trimethylpropan-1-amine"),
]

# 真实 A9 数据（分子仍夹其它未修缺陷：烯炔 (2E)、稠环位次、bis/di、N- 括号等），
# 仅断言本域 N-前缀段正确、无假数字。
nprefix_locant_mixed__SUBSTR_CASES = [
    # must_contain / must_not_contain（在 r.en 中）
    ("CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12.Cl",
     "N,6,6-trimethyl", "0,6,6"),
    ("COc1ccc2c(c1)N(C[C@H](C)CN(C)C)c1ccccc1S2",
     "N,N,2-trimethylpropan-1-amine", "0,0,2"),
    ("Cc1ccccc1-c1cc(N2CCN(C)CC2)ncc1N(C)C(=O)C(C)(C)c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
     "N,2-dimethyl", "1,2-dimethyl"),
    ("Cl.CN(C(=O)C=1SC=CC1C)C1CCNCC1",
     "N,3-dimethyl", "2,3-dimethyl"),
]


@pytest.mark.parametrize("smiles,en", nprefix_locant_mixed__FULL_CASES)
def test_mixed_n_stem_full(smiles: str, en: str) -> None:
    """纯胺/酰胺混合组整串命名：N 在数字前、无假数字。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,must,mustnot", nprefix_locant_mixed__SUBSTR_CASES)
def test_mixed_n_stem_real(smiles: str, must: str, mustnot: str) -> None:
    """真实近错分子（含旁系缺陷）仅校验 N-前缀段与假数字消失。"""
    en = SMILESNNamer().name(smiles).en
    assert must in en
    assert mustnot not in en


# ==========================================================================
# 合并自 test_bridge_suffix_parens.py
# IUPAC: P-14.4,P-16.5.1,P-63.2.2.1,P-65.6.3.2.3
# Layer: L3,L5
#
# O/S/N 桥取代基前缀的括注边界与嵌套层级。
#
# 规则（gold/ChEBI 与 Blue Book 一致）：
# - 前端 R 为自带围栏的复合取代基时，围栏只括前端，-oxy/-sulfanyl/-amino 留在括号外：
#   `...oxan-2-yl]oxymethyl`、`...cyclohex-2-en-1-yl]amino]cyclohexyl`、`...chromen-7-yl]oxy`。
#   依据 P-63.2.2.1.1 的 `(pyridin-2-yl)oxy`、P-65.6.3.2.3 的
#   `3-[(pyridine-3-carbonyl)oxy]propanoic acid`。
# - 前端为简单基（methyl/benzyl/直链 propan-2-yl）或 retained 单词酰基（acetyl/benzoyl/
#   hexadecanoyl）时整段融合平铺：`benzylamino`、`propan-2-yloxy`、`hexadecanoyloxy`。
# - 前端围栏已是方括号且桥为 -amino 时整体再括一层（P-16.5.2 嵌套标记）：
#   `[[(1S)-1-carboxy-3-phenylpropyl]amino]propanoyl`。
# - 中文侧与英文侧同形（括号闭在前端「基」后，氧基/硫基/氨基留括号外）。
#
# 每个期望值均取自 benchmarks/merged_benchmark.json 的同分子 gold（见行末 id）。
# ==========================================================================
bridge_suffix_parens__CASES = [
    # 糖苷 O 桥 + 甲基：前端括号闭在 -yl 后，oxy 与 methyl 融合在外
    (
        "C[C@@H]1O[C@H](OC[C@H]2O[C@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O",
        "(2S,3R,4S,5S,6R)-6-[[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]"
        "oxane-2,3,4,5-tetrol",
        None,
        "chebi-2801",
    ),
    # N 桥复合前端：前端方括号 + amino 整体再括一层（[[X]amino]cyclohexyl）
    (
        "OCC1=C[C@H](N[C@H]2C[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)"
        "[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
        "(2R,3R,4S,5S,6R)-2-[(1R,2R,3S,4S,6R)-2,3-dihydroxy-6-(hydroxymethyl)-4-"
        "[[(1S,4R,5S,6S)-4,5,6-trihydroxy-3-(hydroxymethyl)cyclohex-2-en-1-yl]amino]cyclohexyl]"
        "oxy-6-(hydroxymethyl)oxane-3,4,5-triol",
        None,
        "chebi-2376",
    ),
    # O 桥复合前端（无立体描述符）：括号闭在 -yl 后（P-63.2.2.1.1 的 (pyridin-2-yl)oxy）
    (
        "O1C(=CC=C1)CNC(COC=1C=CC=C2C=CC(=NC12)N1CCCC1)=O",
        "N-(furan-2-ylmethyl)-2-(2-pyrrolidin-1-ylquinolin-8-yl)oxyacetamide",
        "N-(呋喃-2-基甲基)-2-(2-吡咯烷-1-基喹啉-8-基)氧基乙酰胺",
        "tiers-102983",
    ),
    # 苄基型前端融合平铺（甲基直接接桥后缀）
    ("OCCNCc1ccccc1", "2-(benzylamino)ethanol", "2-(苄氨基)乙醇", None),
    # 直链 -yl 前端平铺
    ("BrC=1C=C(C=O)C=C(C1)OC(C)C", "3-bromo-5-propan-2-yloxybenzaldehyde", None, None),
    # 磺酰基桥 + 直链 -yl 前端：英文侧平铺（中文侧仍括注）
    (
        "CC(C)S(=O)(=O)c1ccc(C(N)=O)cc1",
        "4-propan-2-ylsulfonylbenzamide",
        None,
        None,
    ),
]


@pytest.mark.parametrize("smiles,en,zh,gold_id", bridge_suffix_parens__CASES)
def test_bridge_suffix_parens(smiles: str, en: str, zh: str | None, gold_id: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 拆分判据单元测试：(stem, 期望拆分结果)
bridge_suffix_parens__SPLIT_CASES = [
    ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yloxy",
     ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yl", "oxy")),
    ("(8R,9S,13S,14S)-13-methyl-17-oxo-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]"
     "phenanthren-3-yloxy",
     ("(8R,9S,13S,14S)-13-methyl-17-oxo-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]"
      "phenanthren-3-yl", "oxy")),
    ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yloxy", ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yl", "oxy")),
    ("(1S)-1-carboxy-3-phenylpropylamino", ("(1S)-1-carboxy-3-phenylpropyl", "amino")),
    # 不拆：简单保留基、直链 -yl、retained 单词酰基、苄基型前端
    ("benzyloxy", None),
    ("propan-2-yloxy", None),
    ("1,3-dihydroxypropan-2-yloxy", None),
    ("hexadecanoyloxy", None),
    ("(3,4-dichlorophenyl)methylamino", None),
    ("(1-methylimidazol-2-yl)sulfanylmethyl", None),
]


@pytest.mark.parametrize("stem,expected", bridge_suffix_parens__SPLIT_CASES)
def test_split_bridge_suffix(stem: str, expected: tuple[str, str] | None) -> None:
    assert _split_bridge_suffix(stem) == expected


@pytest.mark.parametrize(
    "base,suf,expected",
    [
        ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yl", "oxy", True),
        ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yl", "oxy", True),
        ("(1S)-1-carboxy-3-phenylpropyl", "amino", True),
        ("propan-2-yl", "oxy", False),
        ("1,3-dihydroxypropan-2-yl", "oxy", False),
        ("hexadecanoyl", "oxy", False),
        ("(9Z)-octadec-9-enoyl", "oxy", True),  # 带位次的系统酰基名取 P-65.6.3.2.3 的 [(X)oxy] 式
        ("[4-methyl-2-(trifluoromethyl)phenyl]methyl", "sulfanyl", False),
    ],
)
def test_front_needs_enclosure(base: str, suf: str, expected: bool) -> None:
    assert _front_needs_enclosure(base, suf) is expected


# ==========================================================================
# 合并自 test_simple_alkyl_prefix_paren.py
# IUPAC: P-16.5.1.1
# Layer: L3
#
# Linear n-alkyl as substituent prefix: NO parentheses when simple & locant-free.
#
# Regression for C5+ n-alkyl (pentyl/hexyl/heptyl) that miss the anchored table
# (which only registers up to n-butyl) and fall to the recursive cut→radical-yl
# backend. That backend wrongly marked every non-phenyl/non-alkoxy radical-yl
# paren=True, producing "(pentyl)benzene". Simple unsubstituted alkyl stems get
# no parentheses (P-16.5.1.1: only composite/complex prefixes are enclosed).
#
# Negatives: branched (2,2-dimethylpropyl / 2-methylpropyl) and halo-alkyl
# (3-bromopropyl) stems keep their parentheses — they are composite prefixes.
# ==========================================================================
simple_alkyl_prefix_paren__CASES = [
    # positive: locant-free simple n-alkyl — no parentheses
    ("CCCCCCCc1ccccc1", "heptylbenzene", "庚基苯"),
    ("C1(CCCCC)CCCCC1", "pentylcyclohexane", "戊基环己烷"),
    # negative: composite (locant-carrying / halo) prefixes keep parentheses
    ("CC(C)(C)Cc1ccccc1", "(2,2-dimethylpropyl)benzene", "(2,2-二甲基丙基)苯"),
    ("CC(C)Cc1ccccc1", "(2-methylpropyl)benzene", "(2-甲基丙基)苯"),
    ("BrCCCc1ccccc1", "(3-bromopropyl)benzene", "(3-溴丙基)苯"),
]


@pytest.mark.parametrize("smiles,en,zh", simple_alkyl_prefix_paren__CASES)
def test_simple_alkyl_prefix_paren(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_omit_locant_c2_fg.py
# IUPAC: P-14.3.4.4 / P-66 (acetonitrile 省略例) / P-69 L136 (2-substituted ethanol)
# Layer: L5 (assembler_prefixes._omit_sub_locants)
#
# C2 单取代母体 locant 省略收窄：醇/硫醇/胺 FG 端碳有可取代 H，2- 不可省略。
#
# 原逻辑对所有 C2 单取代母体一刀切省略（cyclopropylethanol、aminoethanol、chloroethanol）。
# 需排除 kind ∈ {alcohol, amine, thiol}；而端碳无 H 的腈/酸等仍可省略（(1H-indol-5-yl)acetonitrile）。
# ==========================================================================
omit_locant_c2_fg__FULL_CASES = [
    ("NCCO", "2-aminoethanol", "2-氨基乙醇"),
    ("ClCCO", "2-chloroethanol", "2-氯乙醇"),
    ("C1(CC1)CCO", "2-cyclopropylethanol", "2-环丙基乙醇"),
]

# 仍允许省略（端碳无 H / 对称）—— 防止过度收窄
omit_locant_c2_fg__KEEP_OMIT_CASES = [
    ("N1C=CC2=CC(=CC=C12)CC#N", "(1H-indol-5-yl)acetonitrile"),
]


@pytest.mark.parametrize("smiles,en,zh", omit_locant_c2_fg__FULL_CASES)
def test_c2_fg_keeps_locant(smiles: str, en: str, zh: str) -> None:
    """C2 醇/胺单取代母体：2- 必须保留。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en", omit_locant_c2_fg__KEEP_OMIT_CASES)
def test_c2_fg_still_omits(smiles: str, en: str) -> None:
    """端碳无可取代 H 的腈母体：2- 省略仍成立。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# ==========================================================================
# 合并自 test_tert_pentyl.py
# IUPAC: P-29.3.2
# Layer: L2,L3,L5
#
# 2-methylbutan-2-yl (tert-pentyl / 1,1-dimethylpropyl) branched alkyl.
#
# IUPAC P-29.3.2 systematic branched alkyl prefix. Topology:
#   parent–C(CH3)(CH3)–CH2–CH3
# Shared via L2 `_side_atoms` / L3 `_BRANCH_CHECKS`. L5 wraps leading-locant
# stems in parentheses. tert-butyl and other retained branches must not regress.
# ==========================================================================
tert_pentyl__CASES = [
    # positive: mono 2-methylbutan-2-yl cycloalkane / arene
    (
        "C1CC(C(C)(CC)C)CCC1",
        "(2-methylbutan-2-yl)cyclohexane",
        "(2-甲基丁-2-基)环己烷",
    ),
    (
        "CC(C)(CC)c1ccccc1",
        "(2-methylbutan-2-yl)benzene",
        "(2-甲基丁-2-基)苯",
    ),
    # negative: near-miss retained / linear must not break
    ("CC(C)(C)C1CCCCC1", "tert-butylcyclohexane", "叔丁基环己烷"),
    ("CC(C)(C)c1ccccc1", "tert-butylbenzene", "叔丁基苯"),
    ("CCC", "propane", "丙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", tert_pentyl__CASES)
def test_tert_pentyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("C1CC(C(C)(CC)C)CCC1", "cyclohexane"),
        ("CCC(C)(C)C1CCCCC1", "cyclohexane"),
        ("CC(C)(CC)c1ccccc1", "3,3-dimethylnonane"),
    ],
)
def test_not_silent_or_chain(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) != normalize_en(bad)


@pytest.mark.parametrize(
    "stem,key",
    [
        ("2-methylbutan-2-yl", "methylbutan-2-yl"),
        ("tert-butyl", "butyl"),
        ("sec-butyl", "butyl"),
    ],
)
def test_alkyl_alpha_key_tert_pentyl(stem: str, key: str) -> None:
    assert alkyl_alpha_key(stem) == key


# ==========================================================================
# 合并自 test_branched_alkyl_mode.py
#
# Verify general vs pin mode name switching for branched alkyls.
# ==========================================================================
branched_alkyl_mode__CASES = [
    # (smiles, general_en_substring, pin_en_substring)
    # Must use a parent (benzene) so branched alkyl is extracted as substituent.
    ("CC(C)(C)c1ccccc1", "tert-butyl", "tert-butyl"),  # PIN level: same in both
]


@pytest.mark.parametrize("smiles,gen_token,pin_token", branched_alkyl_mode__CASES)
def test_mode_switching(smiles, gen_token, pin_token):
    r_gen = SMILESNNamer().name(smiles)
    r_pin = SMILESNNamer(name_mode="pin").name(smiles)
    assert r_gen.success and r_pin.success
    assert gen_token in normalize_en(r_gen.en)
    assert pin_token in normalize_en(r_pin.en)


def test_default_is_general():
    """SMILESNNamer() == SMILESNNamer(name_mode='general')."""
    r_default = SMILESNNamer().name("CC(C)CC")
    r_general = SMILESNNamer(name_mode="general").name("CC(C)CC")
    assert normalize_en(r_default.en) == normalize_en(r_general.en)


def test_isopropyl_pin_uses_propan_2_yl():
    r = SMILESNNamer(name_mode="pin").name("CC(C)c1ccccc1")
    assert r.success
    assert "propan-2-yl" in normalize_en(r.en)
    assert "isopropyl" not in normalize_en(r.en)


def test_tert_butyl_unchanged_in_pin():
    """tert-butyl IS a PIN — no change in pin mode."""
    r_gen = SMILESNNamer().name("CC(C)(C)c1ccccc1")
    r_pin = SMILESNNamer(name_mode="pin").name("CC(C)(C)c1ccccc1")
    assert "tert-butyl" in normalize_en(r_gen.en)
    assert "tert-butyl" in normalize_en(r_pin.en)


# ==========================================================================
# 合并自 test_monoalkyl_chain.py
# IUPAC: P-29.3.1
# Layer: L3,L4,L5
#
# Straight-chain monoalkyl (C1–C4) side-chain prefixes on alkane/alcohol parents.
# ==========================================================================
monoalkyl_chain__CASES = [
    # positive: methyl on alcohol / alkane
    ("CC(C)CO", "2-methylpropan-1-ol", "2-甲基丙-1-醇"),
    ("CC(C)(C)O", "2-methylpropan-2-ol", "2-甲基丙-2-醇"),
    ("CCC(C)(C)O", "2-methylbutan-2-ol", None),
    ("CC(O)CC", "butan-2-ol", "丁-2-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", monoalkyl_chain__CASES)
def test_monoalkyl_chain(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
