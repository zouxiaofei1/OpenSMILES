"""L5 names for simple mono-sulfonate ester parents (P-65.3.2)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, to_benzene


_ALKYL_SO_EN = {1: "methanesulfonate", 2: "ethanesulfonate",
                3: "propanesulfonate", 4: "butanesulfonate"}
_ALKYL_SO_ZH = {1: "甲磺酸", 2: "乙磺酸", 3: "丙磺酸", 4: "丁磺酸"}


def _alkyl_o_pair(n: int) -> tuple[str, str] | None:
    """O-side alkyl: EN methyl…; ZH 甲 / 乙 / 丙 / 正丁."""
    from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh
    en = ester_alkyl_en(n)
    if not en:
        return None
    if n == 4:
        return en, "正丁"
    zh = ester_alkyl_zh(n)
    return (en, zh) if zh else None


def _toluene_retained(en: str, zh: str) -> tuple[str, str]:
    """Gold prefers 对甲苯磺酸 over 4-甲基苯磺酸 for p-tosylates."""
    if en == "4-methylbenzene":
        return en, "对甲苯"
    if en == "2-methylbenzene":
        return en, "邻甲苯"
    if en == "3-methylbenzene":
        return en, "间甲苯"
    return en, zh


def _alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SO_EN.get(n), _ALKYL_SO_ZH.get(n)
    return (en, zh) if en and zh else None


def _aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = _toluene_retained(*to_benzene(*side_en_zh(s)))
    return f"{en}sulfonate", f"{zh}磺酸"


def _s_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") == "aryl":
        return _aryl_stem(parent)
    if s.get("kind") == "alkyl":
        return _alkyl_stem(int(s.get("n") or 0))
    return None


def _join(o: tuple[str, str], stem: tuple[str, str]) -> tuple[str, str]:
    return f"{o[0]} {stem[0]}", f"{stem[1]}{o[1]}酯"


def sulfonate_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonate":
        return None
    o = _alkyl_o_pair(int((parent.get("o_side") or {}).get("n") or 0))
    stem = _s_stem(parent)
    return _join(o, stem) if o and stem else None
