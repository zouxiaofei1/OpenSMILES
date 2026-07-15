"""L5 functional-class names for alkyl isocyanate / isothiocyanate (P-61.9)."""
from __future__ import annotations

from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh


def _alkyl_pair(n: int) -> tuple[str, str] | None:
    en, zh = ester_alkyl_en(n), ester_alkyl_zh(n)
    return (en, zh) if en and zh else None


def _iso_pair(n: int, en_tail: str, zh_head: str) -> tuple[str, str] | None:
    alk = _alkyl_pair(n)
    if alk is None:
        return None
    en_a, zh_a = alk
    return f"{en_a} {en_tail}", f"{zh_head}{zh_a}酯"


def isocyanate_names(n: int) -> tuple[str, str] | None:
    return _iso_pair(n, "isocyanate", "异氰酸")


def isothiocyanate_names(n: int) -> tuple[str, str] | None:
    return _iso_pair(n, "isothiocyanate", "异硫氰酸")


def iso_kind_names(kind: str, n: int) -> tuple[str, str] | None:
    if kind == "isocyanate":
        return isocyanate_names(n)
    if kind == "isothiocyanate":
        return isothiocyanate_names(n)
    return None
