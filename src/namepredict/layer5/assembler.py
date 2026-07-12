from __future__ import annotations

from namepredict.layer5.stems import ALCOHOL_EN, ALCOHOL_ZH, ALKANE_EN, ALKANE_ZH
from namepredict.types import NameResult


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


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alcohol":
        return _alcohol_names(n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", False))
    return _alkane_names(n)


def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)


def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})


def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    kind, n = _parent_n(numbered)
    names = _names_for(kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    return _ok(names[0], names[1], time_ms, source)
