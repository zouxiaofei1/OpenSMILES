"""L5 names for monocyclic cycloalkane + one exocyclic carbonyl FG."""
from __future__ import annotations

from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH

_KINDS = frozenset({"cycloalkanecarboxylic"})
_PLAIN = {
    "cycloalkanecarboxylic": ("carboxylic acid", "甲酸"),
}


def _cyclo_stem(n: int) -> tuple[str, str] | None:
    """cyclohexane / 环己烷 (keep 烷; do not zh_stem)."""
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or n < 3:
        return None
    return f"cyclo{en[:-1]}e", f"环{zh}"


def _plain(n: int, en_suf: str, zh_suf: str) -> tuple[str, str] | None:
    stem = _cyclo_stem(n)
    return (f"{stem[0]}{en_suf}", f"{stem[1]}{zh_suf}") if stem else None


def _by_kind(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    return _plain(n, *_PLAIN[kind])


def cyclo_exo_fg_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Dispatch cycloalkane exocyclic FG parent stems (no ring locant on FG)."""
    return _by_kind(kind, n, numbered) if kind in _KINDS else None
