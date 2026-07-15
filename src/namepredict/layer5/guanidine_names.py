"""L5 assembly for simple guanidine functional parent (P-66.4.1.2.1)."""
from __future__ import annotations


def _aryl_en_zh(mol, parent: int, c_idx: int) -> tuple[str, str]:
    from namepredict.layer2.aryl_sub import _phenyl_at, _phenyl_name
    ph = _phenyl_at(mol, c_idx, parent)
    if ph is None:
        return "phenyl", "苯基"
    en, zh, _ = _phenyl_name(mol, ph, c_idx)
    return en, zh


def _wrap(en: str, zh: str) -> tuple[str, str]:
    if en == "phenyl":
        return en, zh
    return f"({en})", f"({zh})"


def _unsub_names() -> tuple[str, str]:
    return "guanidine", "胍"


def _aryl_label(sub: dict, mol) -> tuple[str, str]:
    c, n_idx = sub.get("aryl_c"), sub.get("n_idx")
    if mol is None or c is None or n_idx is None:
        return "phenyl", "苯基"
    return _aryl_en_zh(mol, n_idx, c)


def _mono_aryl_names(parent: dict) -> tuple[str, str]:
    en, zh = _wrap(*_aryl_label(parent.get("sub") or {}, parent.get("mol")))
    return f"1-{en}guanidine", f"1-{zh}胍"


def _as_label(sub: dict, mol) -> tuple[str, str]:
    """Aryl name from S-attached Ph; parent of ring is S."""
    c, s_idx = sub.get("aryl_c"), sub.get("s_idx")
    if mol is None or c is None or s_idx is None:
        return "phenyl", "苯基"
    return _aryl_en_zh(mol, s_idx, c)


def _arylsulfonyl_names(parent: dict) -> tuple[str, str]:
    en, zh = _as_label(parent.get("sub") or {}, parent.get("mol"))
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
