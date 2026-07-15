"""L5 names for simple arylboronic acids (P-68.1).

Parent: phenylboronic acid / 苯基硼酸 (unsub ZH: 苯硼酸).
Ring prefixes reuse benzene-style build_prefix (halo / methyl / CF3).
"""
from __future__ import annotations


def _subs(numbered: dict) -> list:
    return list(numbered.get("substituents") or [])


def _paren_cf3_en(en: str) -> str:
    """Multi-sub gold uses 3-bromo-5-(trifluoromethyl)..."""
    if "trifluoromethyl" not in en or "(trifluoromethyl)" in en:
        return en
    return en.replace("trifluoromethyl", "(trifluoromethyl)")


def _prefix_pair(numbered: dict) -> tuple[str, str]:
    from namepredict.layer5.assembler import _build_prefix

    subs = _subs(numbered)
    if not subs:
        return "", ""
    en, zh = _build_prefix(subs, 6, "boronic")
    return _paren_cf3_en(en), zh


def _parent_pair(has_prefix: bool) -> tuple[str, str]:
    # Unsub ZH gold/task: 苯硼酸; sub: 苯基硼酸
    if has_prefix:
        return "phenylboronic acid", "苯基硼酸"
    return "phenylboronic acid", "苯硼酸"


def boronic_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "boronic":
        return None
    pre_en, pre_zh = _prefix_pair(numbered)
    en_p, zh_p = _parent_pair(bool(pre_en or pre_zh))
    en = f"{pre_en}{en_p}" if pre_en else en_p
    zh = f"{pre_zh}{zh_p}" if pre_zh else zh_p
    return en, zh
