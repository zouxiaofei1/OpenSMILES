"""L5 assembly for simple carbamate (P-65)."""
from __future__ import annotations

from namepredict.layer5.stems import ESTER_ALKYL_EN, ESTER_ALKYL_ZH


_N_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _alkyl_pair(n: int | None) -> tuple[str, str] | None:
    if n is None:
        return None
    en, zh = ESTER_ALKYL_EN.get(n), ESTER_ALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkoxy_names(parent: dict) -> tuple[str, str] | None:
    if parent.get("alkoxy_en"):
        return parent["alkoxy_en"], parent.get("alkoxy_zh") or parent["alkoxy_en"]
    return _alkyl_pair(parent.get("alkoxy_n"))


def _n_aryl_piece(parent: dict) -> tuple[str, str]:
    en = parent.get("n_aryl_en") or "phenyl"
    zh = parent.get("n_aryl_zh") or "苯基"
    return (f"({en})", f"({zh})") if en != "phenyl" else ("N-phenyl", "N-苯基")


def _n_alkyl_piece(parent: dict) -> tuple[str, str]:
    n = parent.get("n_alkyl_n")
    en, zh = _N_ALKYL_EN.get(n), _N_ALKYL_ZH.get(n)
    return (f"N-{en}", f"N-{zh}") if en else ("", "")


def _mixed_n(parent: dict) -> tuple[str, str]:
    labels = list(parent.get("n_labels") or [])
    if len(labels) != 2:
        return "", ""
    labels = sorted(labels, key=lambda t: 0 if t.get("kind") == "alk" else 1)
    a, b = labels[0], labels[1]
    return (
        f"{a.get('en') or 'alkyl'}({b.get('en') or 'aryl'})",
        f"{a.get('zh') or '烷基'}({b.get('zh') or '芳基'})",
    )


def _n_piece(parent: dict) -> tuple[str, str]:
    if len(list(parent.get("n_c_idxs") or [])) == 2 and parent.get("n_labels"):
        return _mixed_n(parent)
    if parent.get("n_phenyl") or parent.get("n_aryl_c") is not None:
        return _n_aryl_piece(parent)
    if parent.get("n_benzyl"):
        return "N-benzyl", "N-苄基"
    return _n_alkyl_piece(parent)


def _en_join(a_en: str, n_en: str) -> str:
    if n_en.startswith("N-") or n_en:
        return f"{a_en} {n_en}carbamate"
    return f"{a_en} carbamate"


def _zh_join(a_zh: str, n_zh: str) -> str:
    if n_zh.startswith("N-") or n_zh.startswith("(") or n_zh:
        return f"{n_zh}氨基甲酸{a_zh}酯"
    return f"氨基甲酸{a_zh}酯"


def carbamate_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    alk = _alkoxy_names(parent)
    if not alk:
        return None
    a_en, a_zh = alk
    n_en, n_zh = _n_piece(parent)
    return _en_join(a_en, n_en), _zh_join(a_zh, n_zh)
