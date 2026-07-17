"""L5 names for sat_hetero lactones / lactams (P-65.6.3.5.1 / P-66.1.5.1).

oxolan-2-one / 氧杂环戊烷-2-酮 style: base stem with terminal 'e' dropped for EN,
then -{loc}-one / -{loc}-酮. Locant from oriented chain + ketone_c_idx.
"""
from __future__ import annotations

_STEM = {
    "oxolane": ("oxolane", "氧杂环戊烷"),
    "oxane": ("oxane", "氧杂环己烷"),
    "pyrrolidine": ("pyrrolidine", "吡咯烷"),
    "piperidine": ("piperidine", "哌啶"),
}
_ONE_KINDS = frozenset({"oxolanone", "oxanone", "pyrrolidinone", "piperidinone"})


def _one_loc(numbered: dict) -> int | None:
    parent = numbered.get("parent") or {}
    chain = parent.get("chain") or []
    attach = parent.get("ketone_c_idx")
    if attach is None or attach not in chain:
        return None
    return chain.index(attach) + 1


def _en_one_stem(en: str) -> str:
    return en[:-1] if en.endswith("e") else en


def _pair_from_base(base: str, loc: int) -> tuple[str, str] | None:
    stem = _STEM.get(base)
    if stem is None:
        return None
    en_s, zh_s = stem
    return f"{_en_one_stem(en_s)}-{loc}-one", f"{zh_s}-{loc}-酮"


def sat_hetero_one_names(numbered: dict) -> tuple[str, str] | None:
    """oxolan-2-one / 氧杂环戊烷-2-酮 (etc.)."""
    parent = numbered.get("parent") or {}
    kind = parent.get("kind") or ""
    if kind not in _ONE_KINDS:
        return None
    loc = _one_loc(numbered)
    if loc is None:
        return None
    return _pair_from_base(parent.get("base_kind") or "", loc)
