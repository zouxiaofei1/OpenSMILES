# IUPAC: P-23.2.6.1 / P-23.3.1
# Layer: L5
#
# 桥环 von Baeyer 名组装：环数词头 + 描述符 + 烷词干 + 骨架置换 'a' 前缀。
from __future__ import annotations

import pytest

from opensmiles.layer5.bridged_namer import descriptor_str, ring_count_prefix
from opensmiles.layer5.skeleton_replacement import skeleton_replacement_prefix
from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# 每条 = (SMILES, 期望英文名, 期望中文名)。金标旁证见 benchmarks/merged_benchmark.json。
CASES = [
    ("C1CC2CCC1C2", "bicyclo[2.2.1]heptane", "双环[2.2.1]庚烷"),
    ("C12CNCC(CC1)CC2", "3-azabicyclo[3.2.2]nonane", "3-氮杂双环[3.2.2]壬烷"),
    ("C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.13,7]decane", "三环[3.3.1.13,7]癸烷"),
    # 扭曲烷（C10H16）：萘为其非诱导子图，曾丢环丢氢出 decahydronaphthalene（C10H18）
    ("C12C3CC(C(C1)CC3)CC2", "tricyclo[4.4.0.03,8]decane", "三环[4.4.0.03,8]癸烷"),
    # 0 原子主桥（桥头直键）
    ("O=C1C2C(C2CC1)C(=O)OCC", None, None),
]


@pytest.mark.parametrize("smiles,en,zh", [c for c in CASES if c[1]])
def test_bridged_parent_name(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_bridged_alkene_keeps_double_bond():
    """kind=bridged 必须走链引擎，否则不饱和度静默丢失。"""
    r = SMILESNNamer().name("C1=CC2CCC1C2")
    assert r.success, r
    assert "ene" in r.en, r.en


@pytest.mark.parametrize("smiles,en", [
    # P-23.3.1 的 'a' 前缀 + 环烷，供 >10 元杂单环用（金标 5,8,11-triazacyclononadecane 同型）
    ("C1CCCCCCCCCCCN1", "1-azacyclotridecane"),
    ("C1CCCCCCCCCCCO1", "1-oxacyclotridecane"),
])
def test_large_heteromonocycle_uses_a_prefix(smiles, en):
    """>10 元杂单环取骨架置换前缀；≤10 元保留名路径不受影响。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)


def test_small_heteromonocycle_uses_hw_name_not_a_prefix():
    """≤10 元杂单环取 Hantzsch-Widman 名（P-22.2.2），骨架置换 'a' 前缀不得抢入。"""
    r = SMILESNNamer().name("C1CCCCCCN1")
    assert r.success, r
    assert normalize_en(r.en) == normalize_en("azocane")
    assert r.en != "1-azacyclooctane"  # 'a' 前缀形式只用于 >10 元环


def test_ring_count_prefix_uses_multiplicative_forms():
    """环数词头：2→bi/双，3+ 用数量词（曾错用 en_num_term[:-1] 得 tetrcyclo）。"""
    assert ring_count_prefix(2) == ("bi", "双")
    assert ring_count_prefix(3) == ("tri", "三")
    assert ring_count_prefix(4) == ("tetra", "四")
    assert ring_count_prefix(6) == ("hexa", "六")
    assert ring_count_prefix(8) == ("octa", "八")


def test_descriptor_str_flattens_superscripts():
    """描述符：次级桥位次对直接续在长度数字后（判分口径无 ^{}）。"""
    assert descriptor_str((2, 2, 1), ()) == "2.2.1"
    assert descriptor_str((3, 3, 1, 1), ((3, 7),)) == "3.3.1.13,7"
    assert descriptor_str((10, 2, 1, 0, 0), ((1, 9), (3, 8))) == "10.2.1.01,9.03,8"


def test_skeleton_replacement_prefix_order_and_multiplicity():
    """P-23.3.1：组内位次升序，组间按 F>Cl>Br>I>O>S>…>N 引用。"""
    assert skeleton_replacement_prefix({7: [3]}) == ("3-aza", "3-氮杂")
    assert skeleton_replacement_prefix({8: [3, 14]}) == ("3,14-dioxa", "3,14-二氧杂")
    assert skeleton_replacement_prefix({7: [1], 16: [4]}) == ("4-thia-1-aza", "4-硫杂-1-氮杂")


def test_multiplicative_a_elides_only_before_a():
    """数量词尾 'a' 只在后接元素名以 'a' 开头时省略（金标 tetraza 对 tetraoxa）。"""
    assert skeleton_replacement_prefix({7: [1, 2, 3, 4]}) == ("1,2,3,4-tetraza", "1,2,3,4-四氮杂")
    assert skeleton_replacement_prefix({8: [1, 2, 3, 4]}) == ("1,2,3,4-tetraoxa", "1,2,3,4-四氧杂")
    assert skeleton_replacement_prefix({7: [1, 2], 16: [5]}) == ("5-thia-1,2-diaza", "5-硫杂-1,2-二氮杂")


def test_unsupported_element_yields_none():
    """词表外元素的 'a' 前缀不可臆造。"""
    assert skeleton_replacement_prefix({9: [1]}) is None
