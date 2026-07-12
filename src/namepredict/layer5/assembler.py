from __future__ import annotations

from namepredict.layer5.stems import (
    ACID_EN,
    ACID_ZH,
    ALCOHOL_EN,
    ALCOHOL_ZH,
    ALKANE_EN,
    ALKANE_ZH,
)
from namepredict.types import NameResult

MULT_EN = {
    2: "di",
    3: "tri",
    4: "tetra",
    5: "penta",
    6: "hexa",
    7: "hepta",
    8: "octa",
    9: "nona",
    10: "deca",
}
MULT_ZH = {
    2: "二",
    3: "三",
    4: "四",
    5: "五",
    6: "六",
    7: "七",
    8: "八",
    9: "九",
    10: "十",
}


def _fail(meta: dict | None = None) -> NameResult:
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)


def _pair(en_map: dict, zh_map: dict, n: int) -> tuple[str, str] | None:
    en, zh = en_map.get(n), zh_map.get(n)
    return (en, zh) if en and zh else None


def _alkane_names(n: int) -> tuple[str, str] | None:
    return _pair(ALKANE_EN, ALKANE_ZH, n)


def _alcohol_plain(n: int) -> tuple[str, str] | None:
    return _pair(ALCOHOL_EN, ALCOHOL_ZH, n)


def _alcohol_with_locant(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alcohol_plain(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-2]}-{locant}-ol", f"{zh[0]}-{locant}-醇"


def _omit_oh_locant(n: int, oh_locant: int | None, omit: bool) -> bool:
    return omit or oh_locant is None or (oh_locant == 1 and n <= 2)


def _alcohol_names(n: int, oh_locant: int | None, omit: bool) -> tuple[str, str] | None:
    if _omit_oh_locant(n, oh_locant, omit):
        return _alcohol_plain(n)
    return _alcohol_with_locant(n, oh_locant)


def _acid_names(n: int) -> tuple[str, str] | None:
    return _pair(ACID_EN, ACID_ZH, n)


def _ketone_from_alkane(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-1]}-{locant}-one", f"{zh[0]}-{locant}-酮"


def _ketone_names(n: int, locant: int | None) -> tuple[str, str] | None:
    if locant is None:
        return None
    return _ketone_from_alkane(n, locant)


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alcohol":
        return _alcohol_names(n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", False))
    if kind == "acid":
        return _acid_names(n)
    if kind == "ketone":
        return _ketone_names(n, numbered.get("ketone_locant"))
    return _alkane_names(n)


def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)


def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})


def _group_by_stem(substituents: list) -> dict[str, list]:
    groups: dict[str, list] = {}
    for s in substituents:
        groups.setdefault(s.get("en") or "", []).append(s)
    return groups


def _locant_str(subs: list) -> str:
    locs = sorted(int(s["locant"]) for s in subs if "locant" in s)
    return ",".join(str(x) for x in locs)


def _mult_en(n: int) -> str:
    return MULT_EN.get(n, "")


def _mult_zh(n: int) -> str:
    return MULT_ZH.get(n, "")


def _omit_sub_locants(n_carbons: int, substituents: list) -> bool:
    if n_carbons <= 1:
        return True
    return n_carbons == 2 and len(substituents) == 1


def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    mult = _mult_en(len(subs))
    if omit:
        return f"{mult}{stem}"
    return f"{_locant_str(subs)}-{mult}{stem}"


def _prefix_one_zh(zh_stem: str, subs: list, omit: bool) -> str:
    mult = _mult_zh(len(subs))
    if omit:
        return f"{mult}{zh_stem}"
    return f"{_locant_str(subs)}-{mult}{zh_stem}"


def _sorted_stems(groups: dict[str, list]) -> list[str]:
    return sorted(k for k in groups if k)


def _parts_for_stem(stem: str, subs: list, omit: bool) -> tuple[str, str]:
    zh_stem = subs[0].get("zh") or ""
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit)


def _collect_parts(groups: dict[str, list], omit: bool) -> tuple[list[str], list[str]]:
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for stem in _sorted_stems(groups):
        en_p, zh_p = _parts_for_stem(stem, groups[stem], omit)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts


def _build_prefix(substituents: list, n_carbons: int) -> tuple[str, str]:
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents)
    en_parts, zh_parts = _collect_parts(_group_by_stem(substituents), omit)
    return "-".join(en_parts), "-".join(zh_parts)


def _join_name(prefix: str, parent: str) -> str:
    return f"{prefix}{parent}" if prefix else parent


def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    kind, n = _parent_n(numbered)
    names = _names_for(kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    pre_en, pre_zh = _build_prefix(numbered.get("substituents") or [], n)
    en = _join_name(pre_en, names[0])
    zh = _join_name(pre_zh, names[1])
    return _ok(en, zh, time_ms, source)
