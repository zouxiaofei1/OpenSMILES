"""L5 names for simple mono-sulfonyl chloride parents (P-65.3)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, to_benzene


_ALKYL_SC_EN = {1: "methanesulfonyl chloride", 2: "ethanesulfonyl chloride",
                3: "propanesulfonyl chloride", 4: "butanesulfonyl chloride"}
_ALKYL_SC_ZH = {1: "甲磺酰氯", 2: "乙磺酰氯", 3: "丙磺酰氯", 4: "丁磺酰氯"}


def _alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SC_EN.get(n), _ALKYL_SC_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = to_benzene(*side_en_zh(s))
    return f"{en}sulfonyl chloride", f"{zh}磺酰氯"


_MODE = {
    "alkyl": _alkyl_names,
    "aryl": _aryl_stem,
}


def sulfonyl_chloride_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonyl_chloride":
        return None
    fn = _MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None
