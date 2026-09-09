# IUPAC: P-45.2.2
# Layer: L2,L4
"""并列母体候选按 P-45.2.2 前缀取代基位次集合收窄（P-44.1.1 未决时才生效）。"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent, select_parent_tied
from namepredict.layer4.candidate_keys import prefix_locant_set, suffix_locant_set
from namepredict.namer import SMILESNNamer

# P-45.2.2 正例：后缀位次集合并列、前缀位次集合不同，PIN 取集合更小者。
PIN_CASES = [
    # 文档 P-45.2.2 邻近案例：位次集合 '2,4' 低于 '2,5'。
    ("Oc1c(CCc2cc(O)c(Cl)cc2)cc(Cl)cc1",
     "4-chloro-2-[2-(4-chloro-3-hydroxyphenyl)ethyl]phenol"),
    # 肽：位次集合 '1,2,3,4' 低于 '1,2,3,5'（当前实现误取后者）。
    ("CC(C)[C@H](N=C(O)[C@@H](N=C(O)[C@@H](N=C(O)[C@@H](N)CC(=O)O)C(C)C)C(C)C)C(=O)O",
     "(2S)-2-[[(2S)-2-[[(2S)-2-[[(2S)-2-amino-3-carboxypropanoyl]amino]-3-methylbutanoyl]amino]-3-methylbutanoyl]amino]-3-methylbutanoic acid"),
    ("CC[C@H](C)[C@H](N=C(O)[C@H](CCC(=O)O)N=C(O)[C@@H]1CCCN1)C(=O)O",
     "(2S,3S)-2-[[(2S)-4-carboxy-2-[[(2S)-pyrrolidine-2-carbonyl]amino]butanoyl]amino]-3-methylpentanoic acid"),
]

# 近邻负例：后缀位次集合不同（'2,4,5' vs '3,4,5'），P-44.1.1 先决，P-45.2.2 不得越级。
GUARD_CASES = [
    ("OC[C@H]1OC(O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O",
     "(3R,4S,5R,6R)-6-(hydroxymethyl)-3-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxane-2,4,5-triol"),
]


def test_prefix_locant_set_is_ascending_with_duplicates():
    """位次集合升序排列且保留重复（P-14.3.5）。"""
    numbered = {"substituents": [{"locant": 4}, {"locant": 2}, {"locant": 2}]}
    assert prefix_locant_set(numbered) == ((2, ""), (2, ""), (4, ""))


def test_prefix_locant_set_ignores_locantless_prefixes():
    """近邻负例：无位次前缀（N- 取代基等）不入键，缺 substituents 返回空元组。"""
    assert prefix_locant_set({"substituents": [{"locant": None}, {"locant": 3}]}) == ((3, ""),)
    assert prefix_locant_set({}) == ()


def test_suffix_locant_set_uses_principal_attachment_atoms():
    """P-44.1.1 键取 principal 特征基团附着原子的链上位次，升序去重前保留重复。"""
    numbered = {"parent": {"chain": [10, 11, 12], "principal_attachment_atoms": [12, 10]}}
    assert suffix_locant_set(numbered) == ((1, ""), (3, ""))
    assert suffix_locant_set({"parent": {"chain": [1, 2]}}) == ()


def test_select_parent_tied_head_matches_select_parent():
    """并列组首位与 select_parent 一致（裁决无法决出时行为不变）。"""
    info = analyze(preprocess("Oc1c(CCc2cc(O)c(Cl)cc2)cc(Cl)cc1"))
    head = select_parent_tied(info)[0]
    selected = select_parent(info)
    assert head["owned_atoms"] == selected["owned_atoms"]


def test_select_parent_tied_keeps_multiple_candidates():
    """存在并列候选时返回整组而非单个。"""
    info = analyze(preprocess("Oc1c(CCc2cc(O)c(Cl)cc2)cc(Cl)cc1"))
    assert len(select_parent_tied(info)) > 1


@pytest.mark.parametrize("smiles,en", PIN_CASES)
def test_tied_candidates_resolved_by_p45_2_2(smiles, en):
    """P-45.2.2 决定并列候选：前缀位次集合更小者胜出。"""
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,en", GUARD_CASES)
def test_suffix_locant_difference_blocks_p45_2_2(smiles, en):
    """近邻负例：后缀位次集合不同时 P-44.1.1 先决，P-45.2.2 不介入，候选顺序不变。"""
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


def test_p45_2_2_tie_keeps_previous_candidate():
    """近邻负例：两侧位次集合全平（无位次前缀）时 P-45.2.2 不介入，候选顺序不变。"""
    result = SMILESNNamer().name("O=C(COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O")
    assert result.success
    assert result.en == "(3S,4R,5R)-3,4,5-trihydroxy-2-oxo-6-phosphonooxyhexyl dihydrogen phosphate"


def test_single_candidate_unchanged():
    """近邻负例：无并列候选的普通分子行为不变。"""
    result = SMILESNNamer().name("OCCO")
    assert result.success
    assert normalize_en(result.en) == normalize_en("ethane-1,2-diol")
