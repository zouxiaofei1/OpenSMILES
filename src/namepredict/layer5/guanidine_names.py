"""L5 assembly for simple guanidine functional parent (P-66.4.1.2.1)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, wrap_aryl


def _unsub_names() -> tuple[str, str]:
    return "guanidine", "胍"


def _mono_aryl_names(parent: dict) -> tuple[str, str]:
    en, zh = wrap_aryl(*side_en_zh(parent.get("sub") or {}))
    return f"1-{en}guanidine", f"1-{zh}胍"


def _arylsulfonyl_names(parent: dict) -> tuple[str, str]:
    en, zh = side_en_zh(parent.get("sub") or {})
    return (
        f"1-({en}sulfonyl)guanidine",
        f"1-({zh}磺酰基)胍",
    )


_MODE = {
    "unsub": lambda p: _unsub_names(),
    "aryl": _mono_aryl_names,
    "arylsulfonyl": _arylsulfonyl_names,
}


def guanidine_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "guanidine":
        return None
    fn = _MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None
