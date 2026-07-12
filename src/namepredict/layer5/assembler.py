from __future__ import annotations

from namepredict.types import NameResult

# Alkane parent stems by carbon count (C1–C10)
_ALKANE_EN = {
    1: "methane",
    2: "ethane",
    3: "propane",
    4: "butane",
    5: "pentane",
    6: "hexane",
    7: "heptane",
    8: "octane",
    9: "nonane",
    10: "decane",
}
_ALKANE_ZH = {
    1: "甲烷",
    2: "乙烷",
    3: "丙烷",
    4: "丁烷",
    5: "戊烷",
    6: "己烷",
    7: "庚烷",
    8: "辛烷",
    9: "壬烷",
    10: "癸烷",
}

# Alcohol retained stems by carbon count (not SMILES equality)
_ALCOHOL_EN = {
    1: "methanol",
    2: "ethanol",
    3: "propanol",
    4: "butanol",
    5: "pentanol",
    6: "hexanol",
    7: "heptanol",
    8: "octanol",
    9: "nonanol",
    10: "decanol",
}
_ALCOHOL_ZH = {
    1: "甲醇",
    2: "乙醇",
    3: "丙醇",
    4: "丁醇",
    5: "戊醇",
    6: "己醇",
    7: "庚醇",
    8: "辛醇",
    9: "壬醇",
    10: "癸醇",
}


def _fail(meta: dict | None = None) -> NameResult:
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)


def _alkane_names(n: int) -> tuple[str, str] | None:
    en, zh = _ALKANE_EN.get(n), _ALKANE_ZH.get(n)
    return (en, zh) if en and zh else None


def _alcohol_plain(n: int) -> tuple[str, str] | None:
    en, zh = _ALCOHOL_EN.get(n), _ALCOHOL_ZH.get(n)
    return (en, zh) if en and zh else None


def _alcohol_with_locant(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alcohol_plain(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-2]}-{locant}-ol", f"{zh[0]}{locant}醇"


def _alcohol_names(n: int, oh_locant: int | None, omit: bool) -> tuple[str, str] | None:
    if omit or oh_locant is None or (oh_locant == 1 and n <= 2):
        return _alcohol_plain(n)
    return _alcohol_with_locant(n, oh_locant)


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alcohol":
        return _alcohol_names(n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", False))
    return _alkane_names(n)


def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)


def assemble(
    numbered: dict,
    *,
    time_ms: float = 0.0,
    source: str = "iupac",
) -> NameResult:
    kind, n = _parent_n(numbered)
    names = _names_for(kind, n, numbered)
    if not names:
        return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})
    return _ok(names[0], names[1], time_ms, source)
