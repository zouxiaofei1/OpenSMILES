"""Retained / multi-substituted benzene parent and prefix helpers (P-22.1.3)."""
from __future__ import annotations


def _is_toluene(numbered: dict) -> bool:
    subs = numbered.get("substituents") or []
    return len(subs) == 1 and subs[0].get("kind") == "alkyl" and subs[0].get("n_carbons") == 1


def _is_anisole(numbered: dict) -> bool:
    subs = numbered.get("substituents") or []
    return len(subs) == 1 and subs[0].get("kind") == "alkoxy" and subs[0].get("n_carbons") == 1


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
    if _is_anisole(numbered):
        return "anisole", "苯甲醚"
    if _is_xylene(numbered):
        return "xylene", "苯"
    return "benzene", "苯"


def benzene_prefix(numbered: dict, build_prefix) -> tuple[str, str]:
    if _is_toluene(numbered) or _is_anisole(numbered):
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
    "benzonitrile": ("benzonitrile", "苯甲腈"),
    "benzoyl_chloride": ("benzoyl chloride", "苯甲酰氯"),
    "pyridine": ("pyridine", "吡啶"),
    "furan": ("furan", "呋喃"),
    "thiophene": ("thiophene", "噻吩"),
    "pyrrole": ("1H-pyrrole", "吡咯"),
    "imidazole": ("1H-imidazole", "咪唑"),
    "pyrimidine": ("pyrimidine", "嘧啶"),
    "pyrazine": ("pyrazine", "吡嗪"),
    "pyridazine": ("pyridazine", "哒嗪"),
    "naphthalene": ("naphthalene", "萘"),
}


def arene_fg_parent_names(kind: str) -> tuple[str, str] | None:
    return _ARENE_FG.get(kind)


def _ester_alkyl_pair(alkoxy_n: int | None) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ESTER_ALKYL_EN, ESTER_ALKYL_ZH
    if alkoxy_n is None:
        return None
    en, zh = ESTER_ALKYL_EN.get(alkoxy_n), ESTER_ALKYL_ZH.get(alkoxy_n)
    return (en, zh) if en and zh else None


def benzoate_parent_names(numbered: dict, build_prefix) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    alkyl = _ester_alkyl_pair(parent.get("alkoxy_n"))
    if not alkyl:
        return None
    pre_en, pre_zh = build_prefix(numbered.get("substituents") or [], 6, "benzoate")
    en = f"{alkyl[0]} {pre_en}benzoate" if pre_en else f"{alkyl[0]} benzoate"
    return en, f"{pre_zh}苯甲酸{alkyl[1]}酯"


def benzenediol_names(locs: list[int] | None) -> tuple[str, str] | None:
    if not locs or len(locs) != 2:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"benzene-{loc}-diol", f"苯-{loc}-二酚"


def _pyridine_cooh_loc(numbered: dict) -> int | None:
    parent = numbered.get("parent") or {}
    chain = parent.get("chain") or []
    attach = parent.get("ring_attach_idx")
    if attach is None or attach not in chain:
        return None
    return chain.index(attach) + 1


def pyridinecarboxylic_names(numbered: dict) -> tuple[str, str] | None:
    loc = _pyridine_cooh_loc(numbered)
    if loc is None:
        return None
    return f"pyridine-{loc}-carboxylic acid", f"吡啶-{loc}-羧酸"


def _pyridine_fg_loc(numbered: dict, key: str) -> int | None:
    parent = numbered.get("parent") or {}
    chain = parent.get("chain") or []
    attach = parent.get(key)
    if attach is None or attach not in chain:
        return None
    return chain.index(attach) + 1


def pyridinamine_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant") or _pyridine_fg_loc(numbered, "amine_c_idx")
    if loc is None:
        return None
    return f"pyridin-{loc}-amine", f"吡啶-{loc}-胺"


def pyridinol_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("oh_locant") or _pyridine_fg_loc(numbered, "oh_c_idx")
    if loc is None:
        return None
    return f"pyridin-{loc}-ol", f"吡啶-{loc}-醇"

