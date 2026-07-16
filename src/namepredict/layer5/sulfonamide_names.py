"""L5 names for simple mono-sulfonamide parents (P-65.3)."""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, to_benzene, wrap_aryl


_ALKYL_SA_EN = {1: "methanesulfonamide", 2: "ethanesulfonamide",
                3: "propanesulfonamide", 4: "butanesulfonamide"}
_ALKYL_SA_ZH = {1: "甲磺酰胺", 2: "乙磺酰胺", 3: "丙磺酰胺", 4: "丁磺酰胺"}


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
    return f"{en}sulfonamide", f"{zh}磺酰胺"


def _n_aryl_alkyl(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    stem = _alkyl_stem(int(s.get("n") or 0))
    if stem is None or n.get("kind") != "aryl":
        return None
    ar_en, ar_zh = wrap_aryl(*side_en_zh(n))
    return f"N-{ar_en}{stem[0]}", f"N-{ar_zh}{stem[1]}"


def _split_pref(en: str, zh: str) -> tuple[str, str, str, str]:
    """Split '2,4-dichlorobenzene' → ('2,4-dichloro', 'benzene', ...)."""
    if en.endswith("benzene"):
        return en[: -len("benzene")], "benzene", zh[: -len("苯")], "苯"
    return "", en, "", zh


def _with_n_pref(pe: str, pz: str, cy_en: str, cy_zh: str, be: str, bz: str):
    head_en = f"{pe}-N-" if pe else "N-"
    head_zh = f"{pz}-N-" if pz else "N-"
    return f"{head_en}{cy_en}{be}sulfonamide", f"{head_zh}{cy_zh}{bz}磺酰胺"


def _n_cyclo_aryl(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    if s.get("kind") != "aryl" or n.get("kind") != "cyclo":
        return None
    cy_en, cy_zh = side_en_zh(n, ("", ""))
    if not cy_en:
        return None
    pe, be, pz, bz = _split_pref(*to_benzene(*side_en_zh(s)))
    return _with_n_pref(pe, pz, cy_en, cy_zh, be, bz)


_MODE = {
    "alkyl": _alkyl_names,
    "aryl": _aryl_stem,
    "n_aryl_alkyl": _n_aryl_alkyl,
    "n_cyclo_aryl": _n_cyclo_aryl,
}


def sulfonamide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonamide":
        return None
    fn = _MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None
