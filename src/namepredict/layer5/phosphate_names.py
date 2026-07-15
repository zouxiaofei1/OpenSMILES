"""L5 names for simple monoalkyl phosphate and alkylphosphonic acid (P-67)."""
from __future__ import annotations

from namepredict.layer5.stems import ESTER_ALKYL_EN, ESTER_ALKYL_ZH


def _alkyl_pair(n: int | None) -> tuple[str, str] | None:
    if n is None:
        return None
    en, zh = ESTER_ALKYL_EN.get(n), ESTER_ALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _zh_phos_alkyl(n: int, zh: str) -> str:
    """甲/乙; else 丙基…; C11+ 十三烷基 (project gold style)."""
    if n <= 2:
        return zh
    return f"{zh}烷基" if n > 10 else f"{zh}基"


def phosphate_names(numbered: dict) -> tuple[str, str] | None:
    """methyl dihydrogen phosphate / 磷酸甲酯; tridecyl… / 磷酸十三烷基酯."""
    parent = numbered.get("parent") or {}
    n = parent.get("alkoxy_n") or parent.get("n_carbons")
    alkyl = _alkyl_pair(n)
    if not alkyl or n is None:
        return None
    en_a, zh_a = alkyl
    return f"{en_a} dihydrogen phosphate", f"磷酸{_zh_phos_alkyl(n, zh_a)}酯"


def phosphonic_names(n: int) -> tuple[str, str] | None:
    """methylphosphonic acid / 甲基膦酸."""
    alkyl = _alkyl_pair(n)
    if not alkyl:
        return None
    en_a, zh_a = alkyl
    return f"{en_a}phosphonic acid", f"{zh_a}基膦酸"


def p_fg_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "phosphate":
        return phosphate_names(numbered)
    return phosphonic_names(n) if kind == "phosphonic" else None
