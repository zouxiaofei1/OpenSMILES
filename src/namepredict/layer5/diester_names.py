"""L5 names for open-chain symmetric dialkyl alkanedioates (P-65.1.1 / P-65.6)."""
from __future__ import annotations


def _oxalate_stem() -> tuple[str, str]:
    return "oxalate", "草酸"


def _sys_dioate_stem(n: int) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH, zh_stem

    if n < 3:
        return None
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh:
        return None
    return f"{en}dioate", f"{zh_stem(zh)}二酸"


def _diacid_stem(n: int) -> tuple[str, str] | None:
    """EN dioate stem + ZH diacid stem (without 二…酯)."""
    return _oxalate_stem() if n == 2 else _sys_dioate_stem(n)


def _alkyl_pair(parent: dict) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ester_alkoxy_pair

    return ester_alkoxy_pair(parent)


def diester_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """di{alkyl} {alkane}dioate / {二酸}二{烷}酯."""
    parent = numbered.get("parent") or {}
    if not parent.get("symmetric"):
        return None
    alkyl, stem = _alkyl_pair(parent), _diacid_stem(n)
    if not alkyl or not stem:
        return None
    en_a, zh_a = alkyl
    en_s, zh_s = stem
    return f"di{en_a} {en_s}", f"{zh_s}二{zh_a}酯"
