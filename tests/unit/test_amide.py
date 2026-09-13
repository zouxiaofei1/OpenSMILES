# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_amide.py: Simple primary unsubstituted alkanamides (–CONH2).
test_alkenamide.py: Open-chain monounsaturated monoamides (alkenamides).
test_benzamide.py: Retained parent name benzamide (Ph–C(=O)–N; P-66.1.1 / 中文 6.5.8.1).
test_simple_carbamate.py: Simple alkyl N-…carbamate (incl. Boc / tert-butyl carbamate).
test_amido_contraction.py: R-C(=O)-NH-* 残基的 amido 收缩（P-66.1.1.4.3 方法 1, PIN）：
test_demoted_amide_to_amino.py: 链端酰胺被羧酸（主基团）挤掉后，降级为 oxo + amino 前缀表达，
test_amide_tautomer_normalize.py: CheBI 酰胺烯醇互变异构归一化：C(O)=N → C(=O)-NH。
"""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# 合并自 test_mono_amide.py
# IUPAC: P-66.1.1
# Layer: L1,L2,L4,L5
#
# Simple primary unsubstituted alkanamides (–CONH2).
#
# Primary amide: carbonyl C with =O and N whose only carbon neighbor is that C.
# Retained formamide/acetamide; C≥3 systematic …amide / …酰胺.
# ==========================================================================
mono_amide__CASES = [
    # positive: retained + straight-chain mono primary amides
    ("C(=O)N", "formamide", "甲酰胺"),
    ("CCCCCC(=O)N", "hexanamide", "己酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_amide__CASES)
def test_mono_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenamide.py
# IUPAC: P-66.1.1 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated monoamides (alkenamides).
#
# Amide is the principal characteristic group (C(=O)N carbon = locant 1); one
# non-aromatic C=C is expressed as -n-enamide / -n-烯酰胺 with the lower
# double-bond carbon locant and (E)/(Z) when stereo is defined (P-31.1).
# ==========================================================================
alkenamide__CASES = [
    ("NC(=O)/C=C/C", "(2E)-but-2-enamide", "(2E)-丁-2-烯酰胺"),
    ("NC(=O)CC=C", "but-3-enamide", "丁-3-烯酰胺"),
    ("NC(=O)/C=C/c1ccccc1", "(2E)-3-phenylprop-2-enamide", "(2E)-3-苯基丙-2-烯酰胺"),
    (r"CCCCCCCC/C=C\CCCCCCCC(=O)N", "(9Z)-octadec-9-enamide", "(9Z)-十八-9-烯酰胺"),
    ("NC(=O)CCC", "butanamide", "丁酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenamide__CASES)
def test_alkenamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzamide.py
# IUPAC: P-66.1.1 / P-65.1.1.1
# Layer: L2,L3,L4,L5
#
# Retained parent name benzamide (Ph–C(=O)–N; P-66.1.1 / 中文 6.5.8.1).
#
# Single amide, single benzene, direct Ph–C(=O)–N. Ring ≤3 simple prefixes
# (halo / n-alkyl·Me / OH / alkoxy / amino / nitro / CF3). N is H or simple
# N-C1–C4 / N,N-dialkyl / N-phenyl (L3 N-extract). Prevent formamide collapse.
# ==========================================================================
benzamide__CASES = [
    # positive: unsubstituted retained parent
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    # positive: ring simple prefixes (attach = 1)
    ("O=C(N)c1ccc(Cl)cc1", "4-chlorobenzamide", "4-氯苯甲酰胺"),
    ("O=C(N)c1ccc(O)cc1", "4-hydroxybenzamide", "4-羟基苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("O=C(N)c1ccc(N)cc1", "4-aminobenzamide", "4-氨基苯甲酰胺"),
    ("O=C(N)c1ccc([N+](=O)[O-])cc1", "4-nitrobenzamide", "4-硝基苯甲酰胺"),
    ("c1ccccc1C=O", "benzaldehyde", "苯甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", benzamide__CASES)
def test_benzamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzamide_not_formamide_collapse() -> None:
    """Ar–CONH2 must not collapse to phenylformamide / formamide."""
    r = SMILESNNamer().name("c1ccccc1C(=O)N")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzamide"
    assert "formamide" not in en
    assert "phenyl" not in en


# ==========================================================================
# 合并自 test_simple_carbamate.py
# IUPAC: P-65
# Layer: L1,L2,L5
#
# Simple alkyl N-…carbamate (incl. Boc / tert-butyl carbamate).
#
# Carbamate R2N–C(=O)–OR is principal FG (not ester). English follows gold:
#   methyl N-methylcarbamate; tert-butyl (3-methoxyphenyl)carbamate.
# ==========================================================================



# ==========================================================================
# 合并自 test_amido_contraction.py
# IUPAC: P-66.1.1.4.3
# Layer: L3,L5
#
# R-C(=O)-NH-* 残基的 amido 收缩（P-66.1.1.4.3 方法 1, PIN）：
# acetyl/formyl/benzoyl 三词收成 acetamido/formamido/benzamido（gold/ChEBI 全量仅此三词；
# 长链/烯酰/被取代苯甲酰/杂环羰酰 gold 保持 acylamino 方法 2——见 test_acyl_oyl.py、
# test_ring_acyl_carbonyl.py 锁定，本规则不得误伤）。中文 gold 对糖类用「乙酰氨基」与
# 本表不完全对齐，故仅精确断言英文。
#
# 生成点：L5 assembler._mononuclear_radical_names 单取代 azane 分支（现 free_to_yl 出
# acetylamino）；括号由 as_substituent._radical_yl_from_sub composite 判定给 (…amino)。
# ==========================================================================
amido_contraction__POS_FRAG = [
    ("*NC(=O)C", "acetamido"),
    ("*NC=O", "formamido"),
    ("*NC(=O)c1ccccc1", "benzamido"),
]


@pytest.mark.parametrize("smiles,en", amido_contraction__POS_FRAG)
def test_amido_fragment(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 片段级近邻负例：不在收缩表内不得误伤（长链酰方法 2、被取代苯甲酰仍 acylamino）。
amido_contraction__NEG_FRAG = [
    ("*NC(=O)CCC", "butanoylamino"),
    ("*NC(=O)c1ccc(Cl)cc1", "4-chlorobenzoylamino"),
    ("*NC(=O)c1ccco1", "furan-2-carbonylamino"),
]


@pytest.mark.parametrize("smiles,en", amido_contraction__NEG_FRAG)
def test_amido_fragment_negative(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例（含 N,N 双乙酰…等倍增场景落 di- 而非 bis-）。
amido_contraction__POS_WHOLE = [
    ("CC(O)=N[C@@H](Cc1ccccc1)C(=O)O", "(2S)-2-acetamido-3-phenylpropanoic acid"),
    ("CC(O)=Nc1ccc(C(=O)O)c(O)c1", "4-acetamido-2-hydroxybenzoic acid"),
    ("O=C(O)CN=CO", "2-formamidoacetic acid"),
    ("O=C(O)CN=C(O)c1ccccc1", "2-benzamidoacetic acid"),
    ("CC(O)=Nc1c(I)c(N=C(C)O)c(I)c(C(=O)O)c1I",
     "3,5-diacetamido-2,4,6-triiodobenzoic acid"),
]


@pytest.mark.parametrize("smiles,en", amido_contraction__POS_WHOLE)
def test_amido_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# ==========================================================================
# 合并自 test_demoted_amide_to_amino.py
# IUPAC: P-41 / P-65.1.2
# Layer: L1,L3,L5
#
# 链端酰胺被羧酸（主基团）挤掉后，降级为 oxo + amino 前缀表达，
# 且 amino/hydroxy 是简单单核前缀，不得带括号。
#
# 回归：NC(=O)… 端酰胺 + COOH 主链，曾错误输出 5-(amino)-…（递归 claim
# 的 paren 规则误伤）；修复 = L1 _arbitrate_parts 把降级酰胺的伯 N 回收
# 为 amines，走正规 amino 前缀通道。
# ==========================================================================
demoted_amide_to_amino__CASES = [
    # 用户例（酮式）：HOOC-CH(OH)-CH2-CH2-C(=O)NH2
    ("NC(=O)CCC(O)C(=O)O", "5-amino-2-hydroxy-5-oxopentanoic acid", "5-氨基-2-羟基-5-氧代戊酸"),
    # 用户例（烯醇/亚胺酸写法，互变异构后应同名）
    ("N=C(O)CCC(O)C(=O)O", "5-amino-2-hydroxy-5-oxopentanoic acid", "5-氨基-2-羟基-5-氧代戊酸"),
    # 无 α-OH 变体
    ("NC(=O)CCCC(=O)O", "5-amino-5-oxopentanoic acid", None),
]


@pytest.mark.parametrize("smiles,en,zh", demoted_amide_to_amino__CASES)
def test_demoted_amide_n_as_amino_no_paren(smiles: str, en: str, zh: str | None) -> None:
    """端酰胺被羧酸挤掉 → amino 走正规前缀，不带括号。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_amide_tautomer_normalize.py
# Layer: L0(tautomer), L1, L2, L5
#
# CheBI 酰胺烯醇互变异构归一化：C(O)=N → C(=O)-NH。
#
# 非芳香中性 C(OH)=N 写法（亚胺酸/烯醇）应被当作酮式酰胺解析与命名；
# 带电 N / O⁻ 阴离子 / 芳香烯醇（吡啶酮类）不归一化。
# ==========================================================================
def amide_tautomer_normalize___canon(smiles: str) -> str:
    """返回 SMILES 对应结构的规范串（含立体）。"""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)


def test_preprocess_normalizes_enol_to_amide() -> None:
    """preprocess 后烯醇与酮式 canonical 相同，且 L1 检为 amide 而非醇。"""
    mol = preprocess("CC(O)=NC")
    assert amide_tautomer_normalize___canon("CC(=O)NC") == Chem.MolToSmiles(mol)
    info = analyze(mol)
    assert info["amides"]
    assert not (info.get("hydroxyls") or [])


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(O)=NC", "N-methylacetamide"),
        ("CCCCN=CO", "N-butylformamide"),
        ("CC(O)=NCCc1cnc[nH]1", "N-[2-(1H-imidazol-5-yl)ethyl]acetamide"),
    ],
)
def test_enol_named_as_amide(smiles: str, expected: str) -> None:
    """烯醇写法命出的酰胺名与其酮式/参考一致。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(expected)


def test_enol_equals_keto_naming() -> None:
    """同一化合物的烯醇与酮式写法得到相同名称。"""
    n = SMILESNNamer()
    for enol, keto in [
        ("CC(O)=NC", "CC(=O)NC"),
        ("N[C@@H](CCC(O)=N)C(=O)O", "N[C@@H](CCC(=O)N)C(=O)O"),
    ]:
        e = n.name(enol)
        k = n.name(keto)
        assert e.success and k.success
        assert normalize_en(e.en) == normalize_en(k.en)


def test_stereo_preserved_across_normalization() -> None:
    """归一化只翻目标键级，不破坏其它手性中心。"""
    enol = preprocess("CC(O)=N[C@@H]1CCC[C@H]1O")
    keto = preprocess("CC(=O)N[C@@H]1CCC[C@H]1O")
    assert Chem.MolToSmiles(enol) == Chem.MolToSmiles(keto)
    r = SMILESNNamer().name("CC(O)=N[C@@H]1CCC[C@H]1O")
    assert r.success
    assert normalize_en(r.en) == normalize_en("N-[(1R,2R)-2-hydroxycyclopentyl]acetamide")


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)NC",  # 已是酮式，不改
        "Oc1ccccn1",  # 芳香烯醇（2-羟基吡啶→吡啶酮类），不改
        "CN(C)CCCSc1ccccc1[NH+]=C(O)/C=C/c1ccccc1",  # 带电 N 亚胺鎓，不改
    ],
)
def test_guard_no_false_normalization(smiles: str) -> None:
    """护栏：非中性/芳香/已酮式结构 preprocess 前后结构不变。"""
    mol = preprocess(smiles)
    assert Chem.MolToSmiles(mol) == amide_tautomer_normalize___canon(smiles)
