"""L5 names for simple mono-sulfonamide parents (P-65.3)."""
from __future__ import annotations


_ALKYL_SA_EN = {1: "methanesulfonamide", 2: "ethanesulfonamide",
                3: "propanesulfonamide", 4: "butanesulfonamide"}
_ALKYL_SA_ZH = {1: "甲磺酰胺", 2: "乙磺酰胺", 3: "丙磺酰胺", 4: "丁磺酰胺"}


def _aryl_en_zh(mol, c_idx: int, parent: int) -> tuple[str, str]:
    from namepredict.layer2.aryl_sub import _phenyl_at, _phenyl_name
    ph = _phenyl_at(mol, c_idx, parent)
    if ph is None:
        return "phenyl", "苯基"
    en, zh, _ = _phenyl_name(mol, ph, c_idx)
    return en, zh


def _to_benzene(en: str, zh: str) -> tuple[str, str]:
    if en.endswith("phenyl"):
        en = en[: -len("phenyl")] + "benzene"
    if zh.endswith("苯基"):
        zh = zh[: -len("苯基")] + "苯"
    return en, zh


def _wrap_aryl(en: str, zh: str) -> tuple[str, str]:
    if en == "phenyl":
        return en, zh
    return f"({en})", f"({zh})"


def _alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SA_EN.get(n), _ALKYL_SA_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _aryl_stem(parent: dict) -> tuple[str, str] | None:
    s, mol = parent.get("s_side") or {}, parent.get("mol")
    if mol is None or s.get("c") is None:
        return None
    en, zh = _to_benzene(*_aryl_en_zh(mol, s["c"], parent["s_idx"]))
    return f"{en}sulfonamide", f"{zh}磺酰胺"


def _n_aryl_alkyl(parent: dict) -> tuple[str, str] | None:
    s, n, mol = parent.get("s_side") or {}, parent.get("n_side") or {}, parent.get("mol")
    stem = _alkyl_stem(int(s.get("n") or 0))
    if stem is None or mol is None or n.get("c") is None:
        return None
    ar_en, ar_zh = _wrap_aryl(*_aryl_en_zh(mol, n["c"], parent["n_idx"]))
    return f"N-{ar_en}{stem[0]}", f"N-{ar_zh}{stem[1]}"


def _cyclo_label(size: int) -> tuple[str, str] | None:
    from namepredict.layer2.side_cycloalkyl import _cycloalkyl_names
    return _cycloalkyl_names(size)


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
    s, n, mol = parent.get("s_side") or {}, parent.get("n_side") or {}, parent.get("mol")
    if mol is None or s.get("c") is None:
        return None
    cy = _cyclo_label(int(n.get("size") or 0))
    if cy is None:
        return None
    pe, be, pz, bz = _split_pref(*_to_benzene(*_aryl_en_zh(mol, s["c"], parent["s_idx"])))
    return _with_n_pref(pe, pz, cy[0], cy[1], be, bz)


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
