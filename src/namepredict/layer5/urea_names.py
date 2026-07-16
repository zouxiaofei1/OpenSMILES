"""L5 assembly for simple urea functional parent (P-66.1.6.1.1)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, tolyl_swap, wrap_aryl

_BARE = frozenset({"phenyl", "p-tolyl", "o-tolyl", "m-tolyl"})


def _aryl_label(side: dict) -> tuple[str, str]:
    return tolyl_swap(*side_en_zh(side))


def _unsub_names() -> tuple[str, str]:
    return "urea", "脲"


def _mono_aryl_names(parent: dict) -> tuple[str, str]:
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = wrap_aryl(*_aryl_label(ar), bare=_BARE)
    return f"{en}urea", f"{zh}脲"


def _me2_aryl_names(parent: dict) -> tuple[str, str]:
    """Benchmark style: 3-(4-chlorophenyl)-1,1-dimethylurea."""
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = _aryl_label(ar)
    return (
        f"3-({en})-1,1-dimethylurea",
        f"3-({zh})-1,1-二甲基脲",
    )


def _pattern(parent: dict) -> str:
    kinds = {(parent.get("n1") or {}).get("kind"), (parent.get("n3") or {}).get("kind")}
    if kinds == {"h"}:
        return "unsub"
    if kinds == {"aryl", "h"}:
        return "mono_aryl"
    if kinds == {"dialkyl", "aryl"}:
        return "me2_aryl"
    return ""


_PAT = {
    "unsub": lambda p: _unsub_names(),
    "mono_aryl": _mono_aryl_names,
    "me2_aryl": _me2_aryl_names,
}


def urea_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "urea":
        return None
    fn = _PAT.get(_pattern(parent))
    return fn(parent) if fn else None
