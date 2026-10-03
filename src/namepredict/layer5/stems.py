"""烷烃/官能团母体碳数词干表与生成器（C1–C9999，P-14.2.1）。"""

from __future__ import annotations

from namepredict.constants import HS_NUMBER, MULT_EN, MULT_ZH, en_num_term, zh_numeral

# --- C1–C10 保留 / 系统基干（字节兼容） ---
_ALKANE_EN_BASE = {
    1: "methane", 2: "ethane", 3: "propane", 4: "butane", 5: "pentane",
    6: "hexane", 7: "heptane", 8: "octane", 9: "nonane", 10: "decane",
}
_ALKANE_ZH_BASE = {
    1: "甲烷", 2: "乙烷", 3: "丙烷", 4: "丁烷", 5: "戊烷",
    6: "己烷", 7: "庚烷", 8: "辛烷", 9: "壬烷", 10: "癸烷",
}
_ZH_SUFFIXES = ("酰胺", "酰氯", "硫醇", "烷", "醇", "酸", "醛", "腈", "胺", "酮", "烯", "炔")


def zh_num(n: int) -> str | None:
    """中文数字词干：1–10 用天干（甲…癸），11+ 复用 zh_numeral。"""
    if n < 1:
        return None
    return HS_NUMBER[n-1] if n <= 10 else zh_numeral(n)


def zh_stem(zh_full: str) -> str:
    """去除中文全名末端的官能团/母体后缀（十一烷→十一）。"""
    for s in _ZH_SUFFIXES:
        if zh_full.endswith(s) and len(zh_full) > len(s):
            return zh_full[: -len(s)]
    return zh_full


def _en_stem(n: int) -> str | None:
    """不含 'ane' 的烷烃词干（≥11 去尾 'a'：undec/icos）。"""
    if n in _ALKANE_EN_BASE:
        return _ALKANE_EN_BASE[n][:-3]
    t = en_num_term(n)
    return t[:-1] if t else None


def alkane_en(n: int) -> str | None:
    """烷烃英文全名：词干 + 'ane'。"""
    s = _en_stem(n)
    return f"{s}ane" if s else None


def alkane_zh(n: int) -> str | None:
    """烷烃中文全名：查保留表或由 zh_num 生成。"""
    if n in _ALKANE_ZH_BASE:
        return _ALKANE_ZH_BASE[n]
    z = zh_num(n)
    return f"{z}烷" if z else None


def stem_forms(body_en: str, body_zh: str, stem_en: str, stem_zh: str) -> tuple[tuple[str, str], tuple[str, str]]:
    """(完整名, 裸词干) 两形态：裸词干 = 完整名去掉 en 尾 'ane' / zh 尾 '烷'。

    FG 分支要带烷的完整名（后缀去 e 用），链式词干引擎要裸词干（自行拼 ane/烷 或 a…-triene）。
    桥环、螺环与大环杂单环生成式共用此形态约定。
    """
    return ((f"{body_en}{stem_en}", f"{body_zh}{stem_zh}"),
            (f"{body_en}{stem_en[:-3]}", f"{body_zh}{stem_zh[:-1]}"))

def acid_to_anion_en(en: str) -> str:
    """酸/酚/硫醇转阴离子英文名（-oic acid→-oate；-ol→-olate；-thiol→-thiolate）。"""
    if en.endswith("oic acid"):
        return en[:-8] + "oate"
    if en.endswith("ic acid"):
        return en[:-7] + "ate"
    if en.endswith("thiol"):  # 硫醇盐（P-66.1.1.4）
        return en + "ate"
    if en.endswith("ol"):  # 酚盐/醇盐（P-66.1.1.4）
        return en + "ate"
    return en


def join_anion_names(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """若母体为羧酸阴离子，将酸后缀转为 -ate / 酸根。"""
    if not (numbered.get("parent") or {}).get("anion"):
        return en, zh
    return acid_to_anion_en(en), (zh if zh.endswith("根") else f"{zh}根")


def _metal_prefix(metal: str | None, n: int, mult: dict) -> str | None:
    """金属名加数量前缀（n=1 无前缀；n>1 用 MULT 表）。"""
    if not metal or n < 1:
        return None
    if n == 1:
        return metal
    m = mult.get(n)
    return f"{m}{metal}" if m else None


def _metal_en_prefix(salt: dict) -> str | None:
    """英文金属数量前缀（源自 salt 记录）。"""
    return _metal_prefix(salt.get("metal"), salt.get("n_metal") or 0, MULT_EN)


def _metal_zh_suffix(salt: dict) -> str | None:
    """中文金属数量后缀（源自 salt 记录）。"""
    return _metal_prefix(salt.get("metal_zh"), salt.get("n_metal") or 0, MULT_ZH)


def join_metal_salt_names(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """应用金属盐、卤化物盐或氢卤酸盐后缀（P-71.2/P-71.3）。"""
    from namepredict.constants import BIS_EN, BIS_ZH

    salt = numbered.get("salt") or {}
    if salt.get("metal"):
        pref = _metal_en_prefix(salt)
        suf = _metal_zh_suffix(salt)
        n_org = int(salt.get("n_org") or 1)  # 有机阴离子份数：>1 用 bis 括起（calcium bis(...acetate)）
        if n_org > 1 and pref and en.endswith("ate"):
            return f"{pref} {BIS_EN.get(n_org, '')}({en})", f"{BIS_ZH.get(n_org, '')}({zh}){suf or ''}"
        if en.endswith("azanide"):  # P-72.2.2.2：氮负离子母体同样按金属盐前缀表达（sodium …azanide）
            return (f"{pref} {en}" if pref else en,
                    f"{salt.get('metal_zh')}盐{zh}" if suf else zh)
        return (
            f"{pref} {en}" if pref and en.endswith("ate") else en,
            zh[:-1] + suf if suf and zh.endswith("根") else zh,  # "…酸根/酚根"去尾缀金属名（P-71.2）
        )
    acid_en = salt.get("acid_salt")
    if acid_en:
        return f"{en} {acid_en}", f"{zh}{salt.get('acid_salt_zh') or acid_en}"
    halide = salt.get("halide")
    if halide:
        hz = salt.get("halide_zh") or halide
        return f"{en} {halide}", (zh[:-1] + hz if zh.endswith("酸根") else zh + hz)
    return en, zh
