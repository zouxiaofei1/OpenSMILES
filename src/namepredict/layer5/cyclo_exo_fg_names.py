"""L5 names for monocyclic cycloalkane + one exocyclic carbonyl FG."""
from __future__ import annotations

from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH, ester_alkoxy_pair

_KINDS = frozenset({
    "cycloalkanecarboxylic", "cycloalkanecarbaldehyde", "cycloalkanecarbonitrile",
    "cycloalkanecarboxamide", "cycloalkanecarboxylate",
    "cycloalkanecarbonyl_chloride", "cycloalkanecarbonyl_bromide",
})
_PLAIN = {
    "cycloalkanecarboxylic": ("carboxylic acid", "甲酸"),
    "cycloalkanecarbaldehyde": ("carbaldehyde", "甲醛"),
    "cycloalkanecarbonitrile": ("carbonitrile", "甲腈"),
    "cycloalkanecarboxamide": ("carboxamide", "甲酰胺"),
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


def _carboxylate(n: int, numbered: dict) -> tuple[str, str] | None:
    stem, alkyl = _cyclo_stem(n), ester_alkoxy_pair(numbered.get("parent") or {})
    if not stem or not alkyl:
        return None
    return f"{alkyl[0]} {stem[0]}carboxylate", f"{stem[1]}甲酸{alkyl[1]}酯"


def _halide(n: int, kind: str) -> tuple[str, str] | None:
    if kind.endswith("bromide"):
        return _plain(n, "carbonyl bromide", "甲酰溴")
    return _plain(n, "carbonyl chloride", "甲酰氯")


def _by_kind(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind in _PLAIN:
        return _plain(n, *_PLAIN[kind])
    if kind == "cycloalkanecarboxylate":
        return _carboxylate(n, numbered)
    return _halide(n, kind)


def cyclo_exo_fg_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Dispatch cycloalkane exocyclic FG parent stems (no ring locant on FG)."""
    return _by_kind(kind, n, numbered) if kind in _KINDS else None
