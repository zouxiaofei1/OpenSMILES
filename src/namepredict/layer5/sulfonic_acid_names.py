"""L5 names for simple mono free sulfonic acid / sulfonate parents (P-65.3)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, to_benzene


_ALKYL_SA_EN = {1: "methanesulfonic acid", 2: "ethanesulfonic acid",
                3: "propanesulfonic acid", 4: "butanesulfonic acid"}
_ALKYL_SA_ZH = {1: "甲磺酸", 2: "乙磺酸", 3: "丙磺酸", 4: "丁磺酸"}


def _alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SA_EN.get(n), _ALKYL_SA_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = to_benzene(*side_en_zh(s))
    return f"{en}sulfonic acid", f"{zh}磺酸"


_MODE = {
    "alkyl": _alkyl_names,
    "aryl": _aryl_stem,
}


def sulfonic_acid_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonic_acid":
        return None
    fn = _MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None
