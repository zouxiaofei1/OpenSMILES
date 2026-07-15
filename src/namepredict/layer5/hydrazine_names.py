"""L5 assembly for simple hydrazine functional parent (P-68.3.1.2)."""
from __future__ import annotations

from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh


def _unsub_names() -> tuple[str, str]:
    return "hydrazine", "肼"


def _phenyl_names() -> tuple[str, str]:
    return "phenylhydrazine", "苯肼"


def _alkyl_label(n: int) -> tuple[str, str] | None:
    en, zh = ester_alkyl_en(n), ester_alkyl_zh(n)
    return (en, zh) if en and zh else None


def _sym_dialkyl_names(ns: list[int]) -> tuple[str, str] | None:
    if len(ns) != 2 or ns[0] != ns[1]:
        return None
    lab = _alkyl_label(ns[0])
    if lab is None:
        return None
    en_a, zh_a = lab
    return f"1,1-di{en_a}hydrazine", f"1,1-二{zh_a}肼"


def _pattern(parent: dict) -> str:
    n1 = (parent.get("n1") or {}).get("kind")
    n2 = (parent.get("n2") or {}).get("kind")
    kinds = {n1, n2}
    if kinds == {"h"}:
        return "unsub"
    if kinds == {"phenyl", "h"}:
        return "phenyl"
    if kinds == {"dialkyl", "h"}:
        return "dialkyl"
    return ""


def _dialkyl_names(parent: dict) -> tuple[str, str] | None:
    n1 = parent.get("n1") or {}
    n2 = parent.get("n2") or {}
    side = n1 if n1.get("kind") == "dialkyl" else n2
    return _sym_dialkyl_names(list(side.get("ns") or []))


_PAT = {
    "unsub": lambda _p: _unsub_names(),
    "phenyl": lambda _p: _phenyl_names(),
    "dialkyl": _dialkyl_names,
}


def hydrazine_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "hydrazine":
        return None
    fn = _PAT.get(_pattern(parent))
    return fn(parent) if fn else None
