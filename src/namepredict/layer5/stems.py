"""Carbon-count stem tables/generators for alkanes and FG parents (C1–C35+).

C1–C10: retained/table. C11–C19: semi-systematic. C20+: multiplicative compose.
IUPAC P-14.2.1 / P-21 (icos- preferred over eicos-).
"""

from __future__ import annotations

# --- C1–C10 retained / systematic base (byte-compatible) ---
_ALKANE_EN_BASE = {
    1: "methane", 2: "ethane", 3: "propane", 4: "butane", 5: "pentane",
    6: "hexane", 7: "heptane", 8: "octane", 9: "nonane", 10: "decane",
}
_ALKANE_ZH_BASE = {
    1: "甲烷", 2: "乙烷", 3: "丙烷", 4: "丁烷", 5: "戊烷",
    6: "己烷", 7: "庚烷", 8: "辛烷", 9: "壬烷", 10: "癸烷",
}
_SEMI_EN = {
    11: "undec", 12: "dodec", 13: "tridec", 14: "tetradec", 15: "pentadec",
    16: "hexadec", 17: "heptadec", 18: "octadec", 19: "nonadec",
}
_UNITS = {
    1: "hen", 2: "do", 3: "tri", 4: "tetra", 5: "penta",
    6: "hexa", 7: "hepta", 8: "octa", 9: "nona",
}
_TENS = {
    2: "icos", 3: "triacont", 4: "tetracont", 5: "pentacont",
    6: "hexacont", 7: "heptacont", 8: "octacont", 9: "nonacont",
}
_DIGIT_ZH = "零一二三四五六七八九"
_ZH_SUFFIXES = ("酰胺", "酰氯", "硫醇", "烷", "醇", "酸", "醛", "腈", "胺", "酮", "烯", "炔")


def zh_num(n: int) -> str | None:
    """Chinese cardinal for n in 1..99 (十一…三十五); None if out of range."""
    if n < 1 or n > 99:
        return None
    if n < 10:
        return _DIGIT_ZH[n]
    tens, ones = divmod(n, 10)
    head = "十" if tens == 1 else f"{_DIGIT_ZH[tens]}十"
    return head if ones == 0 else f"{head}{_DIGIT_ZH[ones]}"


def zh_stem(zh_full: str) -> str:
    """Strip terminal FG/parent suffix from Chinese full name (十一烷→十一)."""
    for s in _ZH_SUFFIXES:
        if zh_full.endswith(s) and len(zh_full) > len(s):
            return zh_full[: -len(s)]
    return zh_full


def _compose_en_stem(n: int) -> str | None:
    """C20+ multiplicative stem (icos / henicos / hexacos / triacont / …)."""
    tens, ones = divmod(n, 10)
    t = _TENS.get(tens)
    if not t or ones > 9:
        return None
    if ones == 0:
        return t
    # icos + do/tri/… → docos/tricos (elide i); hen keeps icos; 30+ no elide
    base = "cos" if tens == 2 and ones >= 2 else t
    return f"{_UNITS[ones]}{base}"


def _en_stem(n: int) -> str | None:
    """Alkane stem without 'ane' (meth…dec / undec… / icos…)."""
    if n in _ALKANE_EN_BASE:
        return _ALKANE_EN_BASE[n][:-3]
    if n in _SEMI_EN:
        return _SEMI_EN[n]
    return _compose_en_stem(n) if n >= 20 else None


def alkane_en(n: int) -> str | None:
    s = _en_stem(n)
    return f"{s}ane" if s else None


def alkane_zh(n: int) -> str | None:
    if n in _ALKANE_ZH_BASE:
        return _ALKANE_ZH_BASE[n]
    z = zh_num(n)
    return f"{z}烷" if z else None


def alcohol_en(n: int) -> str | None:
    s = _en_stem(n)
    return f"{s}anol" if s else None


def alcohol_zh(n: int) -> str | None:
    if n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}醇"
    z = zh_num(n)
    return f"{z}醇" if z else None


def acid_en(n: int) -> str | None:
    if n == 1:
        return "formic acid"
    if n == 2:
        return "acetic acid"
    s = _en_stem(n)
    return f"{s}anoic acid" if s else None


def acid_zh(n: int) -> str | None:
    if n == 1:
        return "甲酸"
    if n == 2:
        return "乙酸"
    if n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}酸"
    z = zh_num(n)
    return f"{z}酸" if z else None


def aldehyde_en(n: int) -> str | None:
    if n == 1:
        return "formaldehyde"
    if n == 2:
        return "acetaldehyde"
    s = _en_stem(n)
    return f"{s}anal" if s else None


def aldehyde_zh(n: int) -> str | None:
    if n == 1:
        return "甲醛"
    if n == 2:
        return "乙醛"
    if n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}醛"
    z = zh_num(n)
    return f"{z}醛" if z else None


def amide_en(n: int) -> str | None:
    if n == 1:
        return "formamide"
    if n == 2:
        return "acetamide"
    s = _en_stem(n)
    return f"{s}anamide" if s else None


def amide_zh(n: int) -> str | None:
    if n == 1:
        return "甲酰胺"
    if n == 2:
        return "乙酰胺"
    if n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}酰胺"
    z = zh_num(n)
    return f"{z}酰胺" if z else None


def nitrile_en(n: int) -> str | None:
    if n == 1:
        return "formonitrile"
    if n == 2:
        return "acetonitrile"
    s = _en_stem(n)
    return f"{s}anenitrile" if s else None


def nitrile_zh(n: int) -> str | None:
    if n == 1:
        return "甲腈"
    if n == 2:
        return "乙腈"
    if n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}腈"
    z = zh_num(n)
    return f"{z}腈" if z else None


def acyl_chloride_en(n: int) -> str | None:
    if n == 2:
        return "acetyl chloride"
    s = _en_stem(n)
    return f"{s}anoyl chloride" if s and n >= 3 else None


def acyl_chloride_zh(n: int) -> str | None:
    if n == 2:
        return "乙酰氯"
    if 3 <= n <= 10 and n in _ALKANE_ZH_BASE:
        return f"{_ALKANE_ZH_BASE[n][0]}酰氯"
    z = zh_num(n)
    return f"{z}酰氯" if z and n >= 3 else None


def ester_acyl_en(n: int) -> str | None:
    if n == 1:
        return "formate"
    if n == 2:
        return "acetate"
    s = _en_stem(n)
    return f"{s}anoate" if s else None


def ester_alkyl_en(n: int) -> str | None:
    retained = {
        1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl",
        5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl",
        9: "nonyl", 10: "decyl",
    }
    if n in retained:
        return retained[n]
    s = _en_stem(n)
    return f"{s}yl" if s else None


def ester_alkyl_zh(n: int) -> str | None:
    base = {
        1: "甲", 2: "乙", 3: "丙", 4: "丁", 5: "戊", 6: "己", 7: "庚", 8: "辛",
        9: "壬", 10: "癸",
    }
    return base.get(n) or zh_num(n)


def _fill(fn, lo: int = 1, hi: int = 35) -> dict[int, str]:
    out: dict[int, str] = {}
    for n in range(lo, hi + 1):
        v = fn(n)
        if v:
            out[n] = v
    return out


# Public dict API (C1–C35 filled by generators; C20+ never handwritten)
ALKANE_EN = _fill(alkane_en)
ALKANE_ZH = _fill(alkane_zh)
ALCOHOL_EN = _fill(alcohol_en)
ALCOHOL_ZH = _fill(alcohol_zh)
ACID_EN = _fill(acid_en)
ACID_ZH = _fill(acid_zh)
ALDEHYDE_EN = _fill(aldehyde_en)
ALDEHYDE_ZH = _fill(aldehyde_zh)
AMIDE_EN = _fill(amide_en)
AMIDE_ZH = _fill(amide_zh)
NITRILE_EN = _fill(nitrile_en)
NITRILE_ZH = _fill(nitrile_zh)
ESTER_ACYL_EN = _fill(ester_acyl_en)
ESTER_ALKYL_EN = _fill(ester_alkyl_en)
ESTER_ALKYL_ZH = _fill(ester_alkyl_zh)
ACYL_CHLORIDE_EN = _fill(acyl_chloride_en, lo=2)
ACYL_CHLORIDE_ZH = _fill(acyl_chloride_zh, lo=2)

# P-63.2.2 symmetric dialkyl ether retained: dimethyl…dibutyl ether / 二…基醚
ETHER_SYM_EN = {
    1: "dimethyl ether", 2: "diethyl ether",
    3: "dipropyl ether", 4: "dibutyl ether",
}
ETHER_SYM_ZH = {
    1: "二甲基醚", 2: "二乙基醚", 3: "二丙基醚", 4: "二丁基醚",
}
# Asymmetric alkoxy prefix (alkoxyalkane)
ALKOXY_EN = {1: "methoxy", 2: "ethoxy", 3: "propoxy", 4: "butoxy"}
ALKOXY_ZH = {1: "甲氧基", 2: "乙氧基", 3: "丙氧基", 4: "丁氧基"}
# P-63.2.1 dialkyl sulfide: dimethyl…dibutyl sulfide / 二…硫醚
SULFIDE_SYM_EN = {
    1: "dimethyl sulfide", 2: "diethyl sulfide",
    3: "dipropyl sulfide", 4: "dibutyl sulfide",
}
SULFIDE_SYM_ZH = {
    1: "二甲硫醚", 2: "二乙硫醚", 3: "二丙硫醚", 4: "二丁硫醚",
}
# Alkyl radical for asymmetric alkyl alkyl sulfide (alphabetical EN)
SULFIDE_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
SULFIDE_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
