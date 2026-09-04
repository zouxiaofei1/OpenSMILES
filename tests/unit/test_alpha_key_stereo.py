# IUPAC: P-14.5
# Layer: L3,L5
"""alkyl_alpha_key 立体描述符组剥除：前导 (2S,3R,…)- / (E)- 不参与字母序，避免数字泄漏使糖基前缀错误排最前。

Affiliated: P-14.5 字母数字序、L5 取代基前缀排序、L4 链方向排序。
"""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key

# 回归：带立体前缀的糖基词干须剥除立体组后按真实词干字母排序，
# 落在 hydroxy / hydroxymethyl 之后（曾因 '2S,…' 数字泄漏被排到最前）。
_SUGAR = "(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yloxy"
_RING2_STEMS = ["hydroxy", "hydroxymethyl", _SUGAR]


def test_sugar_key_no_digit_leak() -> None:
    """糖基词干的排序键不得以位次数字开头（曾泄漏 '2S,…'）。"""
    key = alkyl_alpha_key(_SUGAR)
    assert not key[:1].isdigit()
    assert key.startswith("trihydroxy")


def test_sugar_sorts_after_simple_prefixes() -> None:
    """P-14.5：同一母体上简单前缀排在立体糖基前缀之前。"""
    order = sorted(_RING2_STEMS, key=alkyl_alpha_key)
    assert order == ["hydroxy", "hydroxymethyl", _SUGAR]


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
_GLY_METHYL = (
    "[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yloxy]methyl]"
)


def test_bracket_wrapped_stem_no_ascii_leak() -> None:
    """方括号包裹的糖基-甲基词干不得以 '[' 开头，须剥到 trihydroxy… 实质词干。"""
    key = alkyl_alpha_key(_GLY_METHYL)
    assert not key[:1] in "[("
    assert key.startswith("trihydroxy")


def test_simple_prefix_sorts_before_bracket_complex() -> None:
    """P-14.5：简单羟基词干排在方括号包裹的复合糖基前缀之前（chebi-431 期望序）。"""
    order = sorted(["hydroxy", _GLY_METHYL], key=alkyl_alpha_key)
    assert order == ["hydroxy", _GLY_METHYL]
