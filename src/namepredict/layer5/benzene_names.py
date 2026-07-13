"""Retained / multi-substituted benzene parent and prefix helpers (P-22.1.3)."""
from __future__ import annotations


def _is_toluene(numbered: dict) -> bool:
    subs = numbered.get("substituents") or []
    return len(subs) == 1 and subs[0].get("kind") == "alkyl" and subs[0].get("n_carbons") == 1


def _is_xylene(numbered: dict) -> bool:
    subs = numbered.get("substituents") or []
    if len(subs) != 2:
        return False
    return all(s.get("kind") == "alkyl" and s.get("n_carbons") == 1 for s in subs)


def _xylene_locants(numbered: dict) -> str:
    locs = sorted(int(s["locant"]) for s in numbered.get("substituents") or [] if "locant" in s)
    return ",".join(str(x) for x in locs)


def benzene_parent_names(numbered: dict) -> tuple[str, str]:
    if _is_toluene(numbered):
        return "toluene", "甲苯"
    if _is_xylene(numbered):
        return "xylene", "苯"
    return "benzene", "苯"


def benzene_prefix(numbered: dict, build_prefix) -> tuple[str, str]:
    if _is_toluene(numbered):
        return "", ""
    if _is_xylene(numbered):
        locs = _xylene_locants(numbered)
        return f"{locs}-", f"{locs}-二甲基"
    return build_prefix(numbered.get("substituents") or [], 6, "benzene")


_ARENE_FG = {
    "phenol": ("phenol", "苯酚"), "aniline": ("aniline", "苯胺"),
    "benzoic": ("benzoic acid", "苯甲酸"),
    "benzaldehyde": ("benzaldehyde", "苯甲醛"),
    "acetophenone": ("acetophenone", "苯乙酮"),
}


def arene_fg_parent_names(kind: str) -> tuple[str, str] | None:
    return _ARENE_FG.get(kind)


def benzenediol_names(locs: list[int] | None) -> tuple[str, str] | None:
    if not locs or len(locs) != 2:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"benzene-{loc}-diol", f"苯-{loc}-二酚"

