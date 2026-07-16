"""L5 functional-class names for organic carbonate (P-65.6)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh


def _alkyl_pair(n: int) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh
    en, zh = ester_alkyl_en(n), ester_alkyl_zh(n)
    return (en, zh) if en and zh else None


def _aryl_zh_gold(en: str, zh: str) -> str:
    """Map common dialkylphenyl ZH to retained 甲苯-style gold forms."""
    if en == "2,4-dimethylphenyl" and zh == "2,4-二甲基苯基":
        return "2,4-二甲苯基"
    if en.endswith("methylphenyl") and zh.endswith("甲基苯基"):
        return zh[: -len("甲基苯基")] + "甲苯基"
    return zh


def _pick_sides(parent: dict) -> tuple[dict, dict]:
    return parent.get("side1") or {}, parent.get("side2") or {}


def _alkyl_of(s1: dict, s2: dict) -> dict:
    return s1 if s1.get("kind") == "alkyl" else s2


def _aryl_of(s1: dict, s2: dict) -> dict:
    return s1 if s1.get("kind") == "aryl" else s2


def _sym_names(parent: dict) -> tuple[str, str] | None:
    s1, s2 = _pick_sides(parent)
    n = int((_alkyl_of(s1, s2) or {}).get("n") or 0)
    lab = _alkyl_pair(n)
    if lab is None:
        return None
    en_a, zh_a = lab
    return f"di{en_a} carbonate", f"碳酸二{zh_a}酯"


def _alkyl_aryl_names(parent: dict) -> tuple[str, str] | None:
    s1, s2 = _pick_sides(parent)
    alk, ar = _alkyl_of(s1, s2), _aryl_of(s1, s2)
    lab = _alkyl_pair(int(alk.get("n") or 0))
    if lab is None or ar.get("kind") != "aryl":
        return None
    a_en, a_zh = lab
    ar_en, ar_zh = side_en_zh(ar)
    ar_zh = _aryl_zh_gold(ar_en, ar_zh)
    return f"{a_en} {ar_en} carbonate", f"{ar_zh}{a_zh}基碳酸酯"


def carbonate_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "carbonate":
        return None
    mode = parent.get("mode")
    if mode == "sym_dialkyl":
        return _sym_names(parent)
    if mode == "alkyl_aryl":
        return _alkyl_aryl_names(parent)
    return None
