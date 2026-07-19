"""L5 names for dialkyl sulfone parents (P-65.3.1.2)."""
from __future__ import annotations

_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
_SULFONE_EN = "sulfone"
_SULFONE_ZH = "砜"


def _alkyl(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_EN.get(n), _ALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _sulfone_symmetric(a1, a2):
    en = f"di{a1[0]} {_SULFONE_EN}"
    zh = f"二{a1[1]}{_SULFONE_ZH}"
    return en, zh


def _sulfone_asymmetric(a1, a2):
    en = f"{a1[0]} {a2[0]} {_SULFONE_EN}"
    zh = f"{a1[1]}{a2[1]}{_SULFONE_ZH}"
    return en, zh


def _names(parent: dict) -> tuple[str, str] | None:
    ns = parent.get("alkyl_ns") or (0, 0)
    n1, n2 = ns[0], ns[1]
    if n1 < 1 or n2 < 1:
        return None
    a1 = _alkyl(n1)
    a2 = _alkyl(n2)
    if a1 is None or a2 is None:
        return None
    return _sulfone_symmetric(a1, a2) if n1 == n2 else _sulfone_asymmetric(a1, a2)


def sulfone_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfone":
        return None
    return _names(parent)
