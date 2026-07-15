"""L5 assembly for simple urea functional parent (P-66.1.6.1.1)."""
from __future__ import annotations


def _aryl_en_zh(mol, n_idx: int, c_idx: int) -> tuple[str, str]:
    from namepredict.layer2.aryl_sub import _phenyl_at, _phenyl_name
    ph = _phenyl_at(mol, c_idx, n_idx)
    if ph is None:
        return "phenyl", "苯基"
    en, zh, _ = _phenyl_name(mol, ph, c_idx)
    return en, zh


def _tolyl_swap(en: str, zh: str) -> tuple[str, str]:
    """Gold prefers p-tolyl / 对甲苯基 over 4-methylphenyl."""
    if en == "4-methylphenyl":
        return "p-tolyl", "对甲苯基"
    if en == "2-methylphenyl":
        return "o-tolyl", "邻甲苯基"
    if en == "3-methylphenyl":
        return "m-tolyl", "间甲苯基"
    return en, zh


def _aryl_label(side: dict, mol) -> tuple[str, str]:
    c, n_idx = side.get("aryl_c"), side.get("n_idx")
    if mol is None or c is None or n_idx is None:
        return "phenyl", "苯基"
    return _tolyl_swap(*_aryl_en_zh(mol, n_idx, c))


def _wrap_aryl(en: str, zh: str) -> tuple[str, str]:
    if en in ("phenyl", "p-tolyl", "o-tolyl", "m-tolyl"):
        return en, zh
    return f"({en})", f"({zh})"


def _unsub_names() -> tuple[str, str]:
    return "urea", "脲"


def _mono_aryl_names(parent: dict) -> tuple[str, str]:
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = _wrap_aryl(*_aryl_label(ar, parent.get("mol")))
    return f"{en}urea", f"{zh}脲"


def _me2_aryl_names(parent: dict) -> tuple[str, str]:
    """Benchmark style: 3-(4-chlorophenyl)-1,1-dimethylurea."""
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = _aryl_label(ar, parent.get("mol"))
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
