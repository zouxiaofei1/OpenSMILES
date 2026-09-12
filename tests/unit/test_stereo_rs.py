# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_rs_stereo.py: R/S CIP stereodescriptor prefixes on parent-chain chiral carbons.
test_rs_stereo_expand.py: R/S coverage expand: ketone/ester/sat-hetero + ester slot + collapse gate.
test_rs_stereo_kinds.py: R/S CIP on expanded parent kinds: amide, nitrile, aldehyde, thiol, diacid.
test_ring_rs_stereo.py: 饱和杂环保留母体(sp3)上的 CIP R/S：pyrrolidine/piperidine/morpholine/piperazine/oxolane/oxane 手性中心。
test_alpha_key_stereo.py: alkyl_alpha_key 立体描述符组剥除：前导 (2S,3R,…)- / (E)- 不参与字母序，避免数字泄漏使糖基前缀错误排最前。
test_p144j_stereo_numbering.py: P-14.4(j)：镜像/反向等价编号平局用 CIP 立体描述符破局——较低位次赋予 R/M/r（而非 S/P/s），
test_stereo_accuracy.py: 立体命名正确比例门禁：复用 benchmarks.stereo_benchmark 抽取 benchmark 超长含立体名称
"""
from __future__ import annotations

import json
import pytest

from benchmarks.stereo_benchmark import run_rows, stereo_tokens
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import alkyl_alpha_key, normalize_en, normalize_zh
from pathlib import Path
from rdkit import Chem

# ==========================================================================
# 合并自 test_rs_stereo.py
# IUPAC: P-92 / P-93
# Layer: L4,L5
#
# R/S CIP stereodescriptor prefixes on parent-chain chiral carbons.
# ==========================================================================
rs_stereo__CASES = [

    (
        "CC[C@H](C)C(=O)[O-]",
        "(2S)-2-methylbutanoate",
        "(2S)-2-甲基丁酸根",
    ),
    (
        "NCC[C@@H](O)C[C@H](N)C(=O)O",
        "(2S,4R)-2,6-diamino-4-hydroxyhexanoic acid",
        "(2S,4R)-2,6-二氨基-4-羟基己酸",
    ),
    (
        "O=C(O)C[C@H](O)CCCCCCCO",
        "(3R)-3,10-dihydroxydecanoic acid",
        "(3R)-3,10-二羟基癸酸",
    ),
    (
        "CCCCCCCCCCCC[C@H](O)C(=O)[O-]",
        "(2S)-2-hydroxytetradecanoate",
        "(2S)-2-羟基十四酸根",
    ),
    (
        "C[C@@H](O)CCCCCCC(=O)O",
        "(8R)-8-hydroxynonanoic acid",
        "(8R)-8-羟基壬酸",
    ),

    # 取代基/自由基母体：手性环成为取代基前缀时也应携带自身 R/S
    # （修复1：implicit-H 的 [C@]/[C@@] 被 RDKit 判为 3 配位而非立体中心）
    (
        "O[C@@]1[C@@](C*)OC(O)[C@](N)[C@]1O",
        "[(2R,3R,4S,5R)-5-amino-3,4,6-trihydroxyoxan-2-yl]methyl",
        "[(2R,3R,4S,5R)-5-氨基-3,4,6-三羟基氧杂环己烷-2-基]甲基",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", rs_stereo__CASES)
def test_rs_stereo(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_rs_stereo_expand.py
# IUPAC: P-92 / P-93
# Layer: L5
#
# R/S coverage expand: ketone/ester/sat-hetero + ester slot + collapse gate.
# ==========================================================================
rs_stereo_expand__CASES = [
    # ester: RS between alkyl and acyl (not before methyl)
    (
        "COC(=O)[C@@H](N)CC(C)C",
        "methyl (2S)-2-amino-4-methylpentanoate",
        "(2S)-2-氨基-4-甲基戊酸甲酯",
    ),
    # ketone multi-center
    (
        "O=C(CO)[C@@H](O)[C@H](O)CO",
        "(3S,4R)-1,3,4,5-tetrahydroxypentan-2-one",
        "(3S,4R)-1,3,4,5-四羟基戊-2-酮",
    ),
    # ketone single center with locant
    (
        "CC(=O)[C@H](O)c1ccccc1",
        "(1R)-1-hydroxy-1-phenylpropan-2-one",
        "(1R)-1-羟基-1-苯基丙-2-酮",
    ),

    # regression E/Z + R/S merge
    (
        "C[C@@H](O)CCCC/C=C/C(=O)O",
        "(2E,8R)-8-hydroxynon-2-enoic acid",
        "(2E,8R)-8-羟基壬-2-烯酸",
    ),
    # negatives: no spurious R/S
    ("CCCCCCCCCCCC(=O)O", "dodecanoic acid", "十二酸"),
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", rs_stereo_expand__CASES)
def test_rs_stereo_expand(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_rs_stereo_kinds.py
# IUPAC: P-92 / P-93
# Layer: L5
#
# R/S CIP on expanded parent kinds: amide, nitrile, aldehyde, thiol, diacid.
# ==========================================================================
rs_stereo_kinds__CASES = [
    # positives — RDKit CIP calibrated
    (
        "CC[C@H](C)C(=O)N",
        "(2S)-2-methylbutanamide",
        "(2S)-2-甲基丁酰胺",
    ),
    (
        "CC[C@H](C)C#N",
        "(2S)-2-methylbutanenitrile",
        "(2S)-2-甲基丁腈",
    ),
    (
        "CC[C@H](C)C=O",
        "(2S)-2-methylbutanal",
        "(2S)-2-甲基丁醛",
    ),
    (
        "CC[C@H](C)S",
        "(2S)-butane-2-thiol",
        "(2S)-丁-2-硫醇",
    ),
    (
        "N[C@H](CCCC(=O)O)C(=O)O",
        "(2R)-2-aminohexanedioic acid",
        "(2R)-2-氨基己二酸",
    ),
    (
        "FC1=CC=C(C=C1)[C@@H](CC(=O)N)C=C",
        "(3S)-3-(4-fluorophenyl)pent-4-enamide",
        "(3S)-3-(4-氟苯基)戊-4-烯酰胺",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", rs_stereo_kinds__CASES)
def test_rs_stereo_kinds(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_ring_rs_stereo.py
# IUPAC: P-92 / P-93
# Layer: L5
#
# 饱和杂环保留母体(sp3)上的 CIP R/S：pyrrolidine/piperidine/morpholine/piperazine/oxolane/oxane 手性中心。
# 环母体 kind 不在链式 FG 白名单(rs_fgs)，曾整体丢掉 R/S 前缀。
# ==========================================================================
ring_rs_stereo__CASES = [
    # piperidine 族：2 位(杂原子邻位)与 3 位手性都覆盖
    (
        "C[C@H]1CCCCN1",
        "(2S)-2-methylpiperidine",
        "(2S)-2-甲基哌啶",
    ),
    (
        "C[C@@H]1CCCCN1",
        "(2R)-2-methylpiperidine",
        "(2R)-2-甲基哌啶",
    ),
    (
        "C[C@H]1CNCCC1",
        "(3R)-3-methylpiperidine",
        "(3R)-3-甲基哌啶",
    ),
    # 其余饱和杂环 scaffold，各取一方向
    (
        "C[C@@H]1CCCN1",
        "(2R)-2-methylpyrrolidine",
        "(2R)-2-甲基吡咯烷",
    ),
    (
        "C[C@H]1CCCO1",
        "(2S)-2-methyloxolane",
        "(2S)-2-甲基四氢呋喃",
    ),
    (
        "C[C@@H]1CCCCO1",
        "(2R)-2-methyloxane",
        "(2R)-2-甲基氧杂环己烷",
    ),
    (
        "C[C@@H]1COCCN1",
        "(3R)-3-methylmorpholine",
        "(3R)-3-甲基吗啉",
    ),
    (
        "C[C@@H]1CNCCN1",
        "(2R)-2-methylpiperazine",
        "(2R)-2-甲基哌嗪",
    ),
    # negatives — 无手性中心不伪造 R/S
    ("C[C@H]1CCCCC1", "methylcyclohexane", "甲基环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", ring_rs_stereo__CASES)
def test_ring_rs_stereo(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alpha_key_stereo.py
# IUPAC: P-14.5
# Layer: L3,L5
#
# alkyl_alpha_key 立体描述符组剥除：前导 (2S,3R,…)- / (E)- 不参与字母序，避免数字泄漏使糖基前缀错误排最前。
#
# Affiliated: P-14.5 字母数字序、L5 取代基前缀排序、L4 链方向排序。
# ==========================================================================
alpha_key_stereo___SUGAR = "(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yloxy"
alpha_key_stereo___RING2_STEMS = ["hydroxy", "hydroxymethyl", alpha_key_stereo___SUGAR]


def test_sugar_key_no_digit_leak() -> None:
    """糖基词干的排序键不得以位次数字开头（曾泄漏 '2S,…'）。"""
    key = alkyl_alpha_key(alpha_key_stereo___SUGAR)
    assert not key[:1].isdigit()
    assert key.startswith("trihydroxy")


def test_sugar_sorts_after_simple_prefixes() -> None:
    """P-14.5：同一母体上简单前缀排在立体糖基前缀之前。"""
    order = sorted(alpha_key_stereo___RING2_STEMS, key=alkyl_alpha_key)
    assert order == ["hydroxy", "hydroxymethyl", alpha_key_stereo___SUGAR]


def test_bare_ez_and_rs_lead_stripped() -> None:
    """裸 E/Z 与带位次 R/S 前导立体组都剥除（含其后的连字符）。"""
    assert alkyl_alpha_key("(E)-but-2-en-1-yl") == "but-2-en-1-yl"
    assert alkyl_alpha_key("(2R)-butan-2-yl") == "butan-2-yl"


def test_non_stereo_leading_parens_unaffected() -> None:
    """非立体形状的前导括号组（如 (methylthio)-）不当作立体组剥除。"""
    assert alkyl_alpha_key("(methylthio)methyl") == "methylthio)methyl"


# 回归：chebi-431（N-乙酰氨基糖，C6 连 O-糖基甲基桥）。复合前缀整体被方括号
# 包裹成 '[[(2R,…)-…oxan-2-yloxy]methyl]'，其前导 '[' 不参与字母序，
# 曾因 ASCII '[' 排于所有字母前而被错误置顶；须剥到实质词干 trihydroxy…。
alpha_key_stereo___GLY_METHYL = (
    "[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yloxy]methyl]"
)


def test_bracket_wrapped_stem_no_ascii_leak() -> None:
    """方括号包裹的糖基-甲基词干不得以 '[' 开头，须剥到 trihydroxy… 实质词干。"""
    key = alkyl_alpha_key(alpha_key_stereo___GLY_METHYL)
    assert not key[:1] in "[("
    assert key.startswith("trihydroxy")


def test_simple_prefix_sorts_before_bracket_complex() -> None:
    """P-14.5：简单羟基词干排在方括号包裹的复合糖基前缀之前（chebi-431 期望序）。"""
    order = sorted(["hydroxy", alpha_key_stereo___GLY_METHYL], key=alkyl_alpha_key)
    assert order == ["hydroxy", alpha_key_stereo___GLY_METHYL]


# ==========================================================================
# 合并自 test_p144j_stereo_numbering.py
# IUPAC: P-14.4(j)
# Layer: L4 numbering
#
# P-14.4(j)：镜像/反向等价编号平局用 CIP 立体描述符破局——较低位次赋予 R/M/r（而非 S/P/s），
# 编号方向不随输入 SMILES 写法漂移。chebi-2118 内消旋 cis-1,7-dimethyl-4-(propan-2-yl)cyclodecane
# 两条镜像编号均给最低位次集 {1,4,7}，须确定性地让 R 中心取 locant 1 → (1R,7S)。
# ==========================================================================
p144j_stereo_numbering___CASES = [
    (
        "CC(C)C1CC[C@H](C)CCC[C@H](C)CC1",  # canonical：R 中心位于环书写起点附近
        "(1R,7S)-1,7-dimethyl-4-propan-2-ylcyclodecane",
    ),
    (
        "[C@H]1(C)CCC[C@@H](C)CCC(C(C)C)CC1",  # 重排原子序：S 中心位于分子开头
        "(1R,7S)-1,7-dimethyl-4-propan-2-ylcyclodecane",
    ),
]


def test_p144j_chebi_2118_meso_ring_rs_locant_stable():
    """同一分子不同输入写法，P-14.4(j) 保证都给出 R 中心在 locant 1 的 (1R,7S)。"""
    namer = SMILESNNamer()
    outputs = [namer.name(s).en for s, _ in p144j_stereo_numbering___CASES]
    assert outputs[0] == outputs[1]  # 不再随原子序漂移
    for s, expected in p144j_stereo_numbering___CASES:
        assert namer.name(s).en == expected


def test_p144j_orient_numbering_direction_invariant():
    """numbering 引擎层面：无论初始 chain 从哪个镜像方向给入，定向后 locant 1 恒为 CIP-R 中心。"""
    smi = "CC(C)C1CC[C@H](C)CCC[C@H](C)CC1"
    mol = preprocess(smi)
    info = analyze(mol)
    parent = select_parent(info)[0]
    parent = finalize_parent_ownership(parent, mol)
    subs = extract_substituents(info, parent)
    base = parent["chain"]
    assert base[0] == 3  # 初始链从异丙基碳起，验证下面的收窄确实改写了起点
    # 环上两个手性中心的 CIP（与输入无关的绝对构型）
    from rdkit.Chem import rdCIPLabeler

    m2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))
    Chem.AssignStereochemistry(m2, force=True, cleanIt=True)
    rdCIPLabeler.AssignCIPLabels(m2)
    cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in m2.GetAtoms() if a.HasProp("_CIPCode")}
    # 多种初始编号方向
    start_from_s = [11, 13, 14, 3, 4, 5, 6, 8, 9, 10]  # S 中心做 locant 1 的镜像方向
    orderings = {
        "base": base,
        "rev": list(reversed(base)),
        "start_from_s": start_from_s,
    }
    for name, chain0 in orderings.items():
        oriented = orient_numbering(
            {"mol": mol, "kind": "alkane", "scaffold_id": "carbocycle", "chain": chain0},
            subs,
        )
        locant1 = oriented[0]
        assert cip[locant1] == "R", f"{name}: locant 1 ({locant1}) 应为 R，实际 {cip[locant1]}"


# ==========================================================================
# 合并自 test_stereo_accuracy.py
# IUPAC: P-92/P-93（立体），benchmark 长名样本
#
# 立体命名正确比例门禁：复用 benchmarks.stereo_benchmark 抽取 benchmark 超长含立体名称
# 样本，只比对立体描述符 token（R/S/E/Z，含 fused 位次如 3aR），忽略整名其它差异（前缀
# 排序、氧桥括号写法等）。用于"即使整名不全对也能量出立体正确比例"，并拦截整体 R/S
# 翻转/丢立体的大回归。可加 -s 查看样本统计。
# ==========================================================================
stereo_accuracy__ROOT = Path(__file__).resolve().parents[2]
stereo_accuracy___DATA = stereo_accuracy__ROOT / "benchmarks" / "merged_benchmark.json"
stereo_accuracy___N_SAMPLE = 20  # 样本量：太长 pytest 慢；约 8s


def stereo_accuracy___sample_rows() -> list[dict]:
    """取英文名含立体 token 且长度最长的一批（长名多为深层糖苷/稠环，立体信息量大）。"""
    data = json.loads(stereo_accuracy___DATA.read_text(encoding="utf-8"))
    rows = [
        r for r in data
        if r.get("english_name") and stereo_tokens(r.get("english_name"))
    ]
    rows.sort(key=lambda r: -len(r.get("english_name") or ""))
    return rows[:stereo_accuracy___N_SAMPLE]


def test_stereo_token_extractor() -> None:
    """token 抽取只认立体描述符，不误吞普通单词。"""
    name = "(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl"
    assert stereo_tokens(name) == ["2S", "3R", "4S", "5S", "6R"]
    assert stereo_tokens("5-[(3aR,4R,5R,6aS)-5-hydroxy-2H-cyclopenta[b]furan]") == [
        "3aR", "4R", "5R", "6aS",
    ]
    assert stereo_tokens("but-3-en-1-ol") == []


@pytest.mark.skipif(not stereo_accuracy___DATA.exists(), reason="benchmark data missing")
def test_stereo_ok_ratio_floor() -> None:
    """长立体名样本上，描述符数一致的分子中"立体集合一致"占比不低于下限。"""
    report = run_rows(stereo_accuracy___sample_rows())  # 复用进程内单例命名器（run-scoped 片段缓存）
    s = report["stats"]
    ratio = s["flip_free"] / (s["count_equal"] or 1)
    print(
        f"stereo sample: rows={s['rows']} ok={s['named_ok']} errored={s['errored']} "
        f"count_equal={s['count_equal']} pure_stereo_ok={s['flip_free']} "
        f"ratio={ratio:.1%}"
    )
    # 下限取当前能力的一个余量，只挡"立体被整体丢弃/翻转"类回归；
    # 立体能力提升时应上调该下限，用 -s 看上述统计。
    assert ratio >= 0.20, (
        f"pure stereo ok ratio {ratio:.1%} below floor; stats={s}"
    )
