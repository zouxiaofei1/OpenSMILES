"""L5 names for simple mono free sulfonic acid / sulfonate parents (P-65.3)."""
from __future__ import annotations


_ALKYL_SA_EN = {1: "methanesulfonic acid", 2: "ethanesulfonic acid",
                3: "propanesulfonic acid", 4: "butanesulfonic acid"}
_ALKYL_SA_ZH = {1: "甲磺酸", 2: "乙磺酸", 3: "丙磺酸", 4: "丁磺酸"}


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
