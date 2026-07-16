"""Shared L5 string helpers for special FG aryl stems (no L2 imports)."""
from __future__ import annotations


def tolyl_swap(en: str, zh: str) -> tuple[str, str]:
    """Gold prefers p-tolyl / 对甲苯基 over 4-methylphenyl."""
    if en == "4-methylphenyl":
        return "p-tolyl", "对甲苯基"
    if en == "2-methylphenyl":
        return "o-tolyl", "邻甲苯基"
    if en == "3-methylphenyl":
        return "m-tolyl", "间甲苯基"
    return en, zh


def to_benzene(en: str, zh: str) -> tuple[str, str]:
    """phenyl → benzene stem swap for arenesulfon* parents."""
    if en.endswith("phenyl"):
        en = en[: -len("phenyl")] + "benzene"
    if zh.endswith("苯基"):
        zh = zh[: -len("苯基")] + "苯"
    return en, zh


def wrap_aryl(en: str, zh: str, bare: frozenset[str] | None = None) -> tuple[str, str]:
    """Parenthesize multi-word aryl unless bare (phenyl / tolyl)."""
    bare = bare or frozenset({"phenyl"})
    if en in bare:
        return en, zh
    return f"({en})", f"({zh})"


def side_en_zh(side: dict, default: tuple[str, str] = ("phenyl", "苯基")) -> tuple[str, str]:
    """Read precomputed en/zh from L2 side dict."""
    en, zh = side.get("en"), side.get("zh")
    if en and zh:
        return en, zh
    return default
