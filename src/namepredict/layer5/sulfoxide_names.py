"""L5 names for open-chain dialkyl sulfoxides (P-63.3)."""
from __future__ import annotations

from namepredict.layer5.stems import SULFIDE_ALKYL_EN, SULFIDE_ALKYL_ZH

# Symmetric: dimethyl sulfoxide / 二甲基亚砜 (gold uses 基)
SULFOXIDE_SYM_EN = {
    1: "dimethyl sulfoxide", 2: "diethyl sulfoxide",
    3: "dipropyl sulfoxide", 4: "dibutyl sulfoxide",
}
SULFOXIDE_SYM_ZH = {
    1: "二甲基亚砜", 2: "二乙基亚砜", 3: "二丙基亚砜", 4: "二丁基亚砜",
}


def _sym_names(n: int) -> tuple[str, str] | None:
    en, zh = SULFOXIDE_SYM_EN.get(n), SULFOXIDE_SYM_ZH.get(n)
    return (en, zh) if en and zh else None


def _asym_names(n1: int, n2: int) -> tuple[str, str] | None:
    en1, en2 = SULFIDE_ALKYL_EN.get(n1), SULFIDE_ALKYL_EN.get(n2)
    zh1, zh2 = SULFIDE_ALKYL_ZH.get(n1), SULFIDE_ALKYL_ZH.get(n2)
    if not en1 or not en2 or not zh1 or not zh2:
        return None
    a, b = sorted([(en1, zh1), (en2, zh2)], key=lambda x: x[0])
    return f"{a[0]} {b[0]} sulfoxide", f"{a[1]}{b[1]}亚砜"


def sulfoxide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    ns = parent.get("alkyl_ns")
    if not ns or len(ns) != 2:
        return None
    n1, n2 = int(ns[0]), int(ns[1])
    if n1 == n2:
        return _sym_names(n1)
    return _asym_names(n1, n2)
