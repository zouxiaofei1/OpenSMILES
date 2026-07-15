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


def _sat_diester(n: int, parent: dict) -> tuple[str, str] | None:
    alkyl, stem = _alkyl_pair(parent), _diacid_stem(n)
    if not alkyl or not stem:
        return None
    en_a, zh_a = alkyl
    en_s, zh_s = stem
    return f"di{en_a} {en_s}", f"{zh_s}二{zh_a}酯"


def _unsat_dioate_stem(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.stereo_ez import _ez_prefix
    from namepredict.layer5.unsat_acid import _unsat_acid_pair

    return _unsat_acid_pair(
        n, numbered.get("ene_locant"), _ez_prefix(numbered),
        "enedioate", "烯二酸", 3,
    )


def _unsat_diester(n: int, numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    alkyl, stem = _alkyl_pair(parent), _unsat_dioate_stem(n, numbered)
    if not alkyl or not stem:
        return None
    en_a, zh_a = alkyl
    en_s, zh_s = stem
    return f"di{en_a} {en_s}", f"{zh_s}二{zh_a}酯"


def diester_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """di{alkyl} {alkane}dioate / {二酸}二{烷}酯; unsat adds ene locant + E/Z."""
    parent = numbered.get("parent") or {}
    if not parent.get("symmetric"):
        return None
    if parent.get("double_bond"):
        return _unsat_diester(n, numbered)
    return _sat_diester(n, parent)
