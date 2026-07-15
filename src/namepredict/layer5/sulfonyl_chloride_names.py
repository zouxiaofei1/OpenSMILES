"""L5 names for simple mono-sulfonyl chloride parents (P-65.3)."""
from __future__ import annotations


_ALKYL_SC_EN = {1: "methanesulfonyl chloride", 2: "ethanesulfonyl chloride",
                3: "propanesulfonyl chloride", 4: "butanesulfonyl chloride"}
_ALKYL_SC_ZH = {1: "甲磺酰氯", 2: "乙磺酰氯", 3: "丙磺酰氯", 4: "丁磺酰氯"}


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
    en, zh = _ALKYL_SC_EN.get(n), _ALKYL_SC_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _aryl_stem(parent: dict) -> tuple[str, str] | None:
    s, mol = parent.get("s_side") or {}, parent.get("mol")
    if mol is None or s.get("c") is None:
        return None
    en, zh = _to_benzene(*_aryl_en_zh(mol, s["c"], parent["s_idx"]))
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
