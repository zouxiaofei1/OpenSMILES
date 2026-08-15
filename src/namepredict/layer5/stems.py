"""烷烃与官能团母体的碳数词干表/生成器（C1–C35+）；C1–C10 保留、C11–C19 半系统、C20+ 倍增组词；IUPAC P-14.2.1 / P-21（icos- 优先于 eicos-）。"""

from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH

# --- C1–C10 保留 / 系统基干（字节兼容） ---
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
_HS_NUMBER = "甲乙丙丁戊己庚辛壬癸"
_DIGIT_ZH = "零一二三四五六七八九"
_ZH_SUFFIXES = ("酰胺", "酰氯", "硫醇", "烷", "醇", "酸", "醛", "腈", "胺", "酮", "烯", "炔")


def zh_num(n: int) -> str | None:
    """Chinese stem digit for n: 1..10 天干(甲…癸); 11..99 数字组合(十一…九十九)."""
    if n < 1:
        return None
    if n <= 10:
        return _HS_NUMBER[n-1]
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
    """C20+ 倍增词干（icos / henicos / hexacos / triacont / …）。"""
    tens, ones = divmod(n, 10)
    t = _TENS.get(tens)
    if not t or ones > 9:
        return None
    if ones == 0:
        return t
    # icos + do/tri/… → docos/tricos（省略 i）；hen 保留 icos；30+ 不省略
    base = "cos" if tens == 2 and ones >= 2 else t
    return f"{_UNITS[ones]}{base}"


def _en_stem(n: int) -> str | None:
    """不含 'ane' 的烷烃词干（meth…dec / undec… / icos…）。"""
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

def acid_to_anion_en(en: str) -> str:
    """酸转阴离子英文名：dodecanoic acid → dodecanoate；acetic acid → acetate。"""
    if en.endswith("oic acid"):
        return en[:-8] + "oate"
    if en.endswith("ic acid"):
        return en[:-7] + "ate"
    return en


def acid_to_anion_zh(zh: str) -> str:
    return zh if zh.endswith("根") else f"{zh}根"


def maybe_anion_names(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """If parent is carboxylate anion, convert acid suffix to -ate / 酸根."""
    if not (numbered.get("parent") or {}).get("anion"):
        return en, zh
    return acid_to_anion_en(en), acid_to_anion_zh(zh)


def _metal_prefix(metal: str | None, n: int, mult: dict) -> str | None:
    if not metal or n < 1:
        return None
    if n == 1:
        return metal
    m = mult.get(n)
    return f"{m}{metal}" if m else None


def _metal_en_prefix(salt: dict) -> str | None:
    return _metal_prefix(salt.get("metal"), salt.get("n_metal") or 0, MULT_EN)


def _metal_zh_suffix(salt: dict) -> str | None:
    return _metal_prefix(salt.get("metal_zh"), salt.get("n_metal") or 0, MULT_ZH)


def _salt_en(en: str, salt: dict) -> str:
    pref = _metal_en_prefix(salt)
    return f"{pref} {en}" if pref and en.endswith("ate") else en


def _salt_zh(zh: str, salt: dict) -> str:
    suf = _metal_zh_suffix(salt)
    if not suf or not zh.endswith("酸根"):
        return zh
    return zh[:-1] + suf


def _acid_salt_suffix(salt: dict) -> tuple[str, str] | None:
    en_s, zh_s = salt.get("acid_salt"), salt.get("acid_salt_zh")
    return (en_s, zh_s or en_s) if en_s else None


def _with_acid_salt(en: str, zh: str, salt: dict) -> tuple[str, str]:
    suf = _acid_salt_suffix(salt)
    if suf is None:
        return en, zh
    return f"{en};{suf[0]}", f"{zh};{suf[1]}"


def maybe_metal_salt_names(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """应用碱金属盐或酸式盐（HCl）后缀。"""
    salt = numbered.get("salt") or {}
    if salt.get("metal"):
        return _salt_en(en, salt), _salt_zh(zh, salt)
    return _with_acid_salt(en, zh, salt)



def _fill(fn, lo: int = 1, hi: int = 35) -> dict[int, str]:
    out: dict[int, str] = {}
    for n in range(lo, hi + 1):
        v = fn(n)
        if v:
            out[n] = v
    return out


# 公共字典 API（C1–C35 由生成器填充；C20+ 从不手写）
ALKANE_EN = _fill(alkane_en)
ALKANE_ZH = _fill(alkane_zh)
