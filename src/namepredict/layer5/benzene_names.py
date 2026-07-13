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
    "pyrazole": ("1H-pyrazole", "吡唑"),
    "pyrimidine": ("pyrimidine", "嘧啶"),
    "pyrazine": ("pyrazine", "吡嗪"),
    "pyridazine": ("pyridazine", "哒嗪"),
    "naphthalene": ("naphthalene", "萘"),
    "indole": ("1H-indole", "吲哚"),
    "aziridine": ("aziridine", "氮杂环丙烷"),
    "oxirane": ("oxirane", "环氧乙烷"),
    "oxolane": ("oxolane", "氧杂环戊烷"),
    "oxane": ("oxane", "氧杂环己烷"),
    "pyrrolidine": ("pyrrolidine", "吡咯烷"),
    "piperidine": ("piperidine", "哌啶"),
    "morpholine": ("morpholine", "吗啉"),
    "piperazine": ("piperazine", "哌嗪"),
    "dioxolane": ("1,3-dioxolane", "1,3-二氧戊环"),
    "dioxane": ("1,4-dioxane", "1,4-二氧六环"),
    "thiolane": ("thiolane", "硫杂环戊烷"),
}


def arene_fg_parent_names(kind: str) -> tuple[str, str] | None:
    return _ARENE_FG.get(kind)


def join_parent_name(prefix: str, parent: str) -> str:
    if not prefix:
        return parent
    return f"{prefix}-{parent}" if parent.startswith("1H-") else f"{prefix}{parent}"


def zh_1h_parent(en_parent: str, zh_parent: str, prefix: str) -> str:
    """Prefix Chinese retained 1H-parents with 1H- when ring is substituted."""
    if not prefix or not en_parent.startswith("1H-") or zh_parent.startswith("1H-"):
        return zh_parent
    return f"1H-{zh_parent}"


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
    return f"pyridine-{loc}-carboxylic acid", f"吡啶-{loc}-甲酸"


# base_kind → (en stem, zh stem, needs_1h)
_H5COOH_STEM = {
    "furan": ("furan", "呋喃", False),
    "thiophene": ("thiophene", "噻吩", False),
    "pyrrole": ("pyrrole", "吡咯", True),
    "imidazole": ("imidazole", "咪唑", True),
    "pyrazole": ("pyrazole", "吡唑", True),
}


def _h5cooh_stem(parent: dict) -> tuple[str, str] | None:
    base = parent.get("base_kind") or ""
    got = _H5COOH_STEM.get(base)
    if got is None:
        return None
    en0, zh0, one_h = got
    if one_h:
        return f"1H-{en0}", f"1H-{zh0}"
    return en0, zh0


def hetero5carboxylic_names(numbered: dict) -> tuple[str, str] | None:
    """furan/thiophene/1H-pyrrole/imidazole/pyrazole-n-carboxylic acid / …-甲酸."""
    parent = numbered.get("parent") or {}
    loc = _pyridine_cooh_loc(numbered)
    stem = _h5cooh_stem(parent)
    if loc is None or stem is None:
        return None
    en_s, zh_s = stem
    return f"{en_s}-{loc}-carboxylic acid", f"{zh_s}-{loc}-甲酸"


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


def _pyrimidinamine_en(pre_en: str, loc: int) -> str:
    return f"{pre_en}pyrimidin-{loc}-amine" if pre_en else f"pyrimidin-{loc}-amine"


def _pyrimidinamine_zh(pre_zh: str, loc: int) -> str:
    if not pre_zh:
        return f"嘧啶-{loc}-胺"
    return f"{loc}-氨基-{pre_zh}嘧啶"


def pyrimidinamine_names(numbered: dict, build_prefix) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant") or _pyridine_fg_loc(numbered, "amine_c_idx")
    if loc is None:
        return None
    pre_en, pre_zh = build_prefix(numbered.get("substituents") or [], 6, "pyrimidinamine")
    return _pyrimidinamine_en(pre_en, loc), _pyrimidinamine_zh(pre_zh, loc)


def pyridinol_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("oh_locant") or _pyridine_fg_loc(numbered, "oh_c_idx")
    if loc is None:
        return None
    return f"pyridin-{loc}-ol", f"吡啶-{loc}-醇"


def pyridine_kind_names(kind: str, numbered: dict, build_prefix) -> tuple[str, str] | None:
    if kind == "pyridinecarboxylic": return pyridinecarboxylic_names(numbered)
    if kind == "pyridinamine": return pyridinamine_names(numbered)
    if kind == "pyrimidinamine": return pyrimidinamine_names(numbered, build_prefix)
    if kind == "pyridinol": return pyridinol_names(numbered)
    return None

