"""Retained / multi-substituted benzene parent and prefix helpers (P-22.1.3)."""
from __future__ import annotations

from namepredict.layer2.kind_registry import parent_names as _registry_parent_names
from namepredict.layer3.substituent_extractor import alkyl_alpha_key


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


def _methoxy_subs(subs: list) -> list:
    return [s for s in subs if s.get("kind") == "alkoxy" and s.get("n_carbons") == 1]


def _is_poly_anisole(numbered: dict) -> bool:
    """≥2 ring subs with exactly one methoxy and no alkyl → zh 苯甲醚 parent."""
    subs = numbered.get("substituents") or []
    if len(subs) < 2 or len(_methoxy_subs(subs)) != 1:
        return False
    return not any(s.get("kind") == "alkyl" for s in subs)


def _anisole_reloc(loc: int, meo: int, rev: bool) -> int:
    d = (meo - loc) % 6 if rev else (loc - meo) % 6
    return d + 1


def _anisole_dir_key(others: list, meo: int, rev: bool) -> tuple:
    """Locants in alphabetical citation order (not sorted-set)."""
    ordered = sorted(
        (s for s in others if "locant" in s),
        key=lambda s: alkyl_alpha_key(s.get("en") or ""),
    )
    return tuple(_anisole_reloc(int(s["locant"]), meo, rev) for s in ordered)


def _anisole_pick_rev(others: list, meo: int) -> bool:
    return _anisole_dir_key(others, meo, True) < _anisole_dir_key(others, meo, False)


def _anisole_renum(others: list, meo_loc: int, rev: bool) -> list:
    return [{**s, "locant": _anisole_reloc(int(s["locant"]), meo_loc, rev)} for s in others]


def _anisole_zh_prefix(numbered: dict, build_prefix) -> str:
    subs = numbered.get("substituents") or []
    meo = _methoxy_subs(subs)[0]
    meo_loc = int(meo["locant"])
    others = [s for s in subs if s is not meo]
    renum = _anisole_renum(others, meo_loc, _anisole_pick_rev(others, meo_loc))
    return build_prefix(renum, 6, "benzene")[1]


def _retained_benzene(numbered: dict) -> tuple[str, str] | None:
    if _is_toluene(numbered):
        return "toluene", "甲苯"
    if _is_anisole(numbered):
        return "anisole", "甲氧基苯"
    if _is_xylene(numbered):
        return "xylene", "苯"
    return ("benzene", "苯甲醚") if _is_poly_anisole(numbered) else None


def benzene_parent_names(numbered: dict) -> tuple[str, str]:
    from namepredict.layer5.iso_arene_names import is_fused_iso, iso_fused_parent
    top = _retained_benzene(numbered)
    if top is not None:
        return top
    return iso_fused_parent(numbered) if is_fused_iso(numbered) else ("benzene", "苯")


def _xylene_prefix(numbered: dict) -> tuple[str, str]:
    locs = _xylene_locants(numbered)
    return f"{locs}-", f"{locs}-二甲基"


def benzene_prefix(numbered: dict, build_prefix) -> tuple[str, str]:
    from namepredict.layer5.iso_arene_names import is_fused_iso, iso_fused_prefix
    if _is_toluene(numbered) or _is_anisole(numbered):
        return "", ""
    if _is_xylene(numbered):
        return _xylene_prefix(numbered)
    if is_fused_iso(numbered):
        return iso_fused_prefix(numbered, build_prefix)
    en_pre, zh_pre = build_prefix(numbered.get("substituents") or [], 6, "benzene")
    return (en_pre, _anisole_zh_prefix(numbered, build_prefix)) if _is_poly_anisole(numbered) else (en_pre, zh_pre)


def arene_fg_parent_names(kind: str) -> tuple[str, str] | None:
    return _registry_parent_names(kind)


def _stereo_lead(parent: str) -> tuple[str, str]:
    """Split leading stereo '(E)-' / '(9Z,11E)-' from parent stem."""
    if not parent.startswith("("):
        return "", parent
    close = parent.find(")-")
    if close < 0:
        return "", parent
    return parent[: close + 2], parent[close + 2 :]


def join_parent_name(prefix: str, parent: str) -> str:
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if stem[:1].isdigit() or stem.startswith("1H-") else f"{prefix}{stem}"
    return f"{stereo}{body}"


def _join_alkyl_stereo(alkyl: str, pre: str, rest: str) -> str:
    """alkyl + stereo-lead rest with substituent prefix after stereo."""
    st, body = _stereo_lead(rest)
    return f"{alkyl} {st}{pre}{body}" if alkyl else f"{st}{pre}{body}"


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str]) -> tuple[str, str]:
    """methyl butanoate + 2-oxo → methyl 2-oxobutanoate; keeps (E)- after alkyl."""
    en, zh = names
    if not pre_en:
        return en, zh
    parts = en.split(" ", 1)
    en = _join_alkyl_stereo(parts[0] if len(parts) == 2 else "", pre_en, parts[-1])
    stz, bodyz = _stereo_lead(zh)
    return en, f"{stz}{pre_zh}{bodyz}" if pre_zh else zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
) -> tuple[str, str]:
    if kind in ("ester", "diester"):
        return join_ester_name(pre[0], pre[1], names)
    en = join_parent_name(pre[0], names[0])
    zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    return en, zh


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
    from namepredict.layer5.stems import ester_alkoxy_pair
    parent = numbered.get("parent") or {}
    alkyl = ester_alkoxy_pair(parent)
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


def benzenediamine_names(locs: list[int] | None) -> tuple[str, str] | None:
    if not locs or len(locs) != 2:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"benzene-{loc}-diamine", f"苯-{loc}-二胺"


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


def pyridinecarbonitrile_names(numbered: dict) -> tuple[str, str] | None:
    loc = _pyridine_cooh_loc(numbered)
    if loc is None:
        return None
    return f"pyridine-{loc}-carbonitrile", f"吡啶-{loc}-甲腈"


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


# base sat_hetero kind → (en stem, zh stem) from _ARENE_STEM table
_SHCOOH_STEM = {
    "pyrrolidine": ("pyrrolidine", "吡咯烷"),
    "piperidine": ("piperidine", "哌啶"),
    "piperazine": ("piperazine", "哌嗪"),
    "morpholine": ("morpholine", "吗啉"),
    "oxolane": ("oxolane", "氧杂环戊烷"),
    "oxane": ("oxane", "氧杂环己烷"),
    "thiolane": ("thiolane", "硫杂环戊烷"),
    "aziridine": ("aziridine", "氮杂环丙烷"),
}


def sat_hetero_carboxylic_names(numbered: dict) -> tuple[str, str] | None:
    """piperidine-n-carboxylic acid / 哌啶-n-甲酸 (etc.)."""
    parent = numbered.get("parent") or {}
    loc = _pyridine_cooh_loc(numbered)
    stem = _SHCOOH_STEM.get(parent.get("base_kind") or "")
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


def benzofuranamine_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant")
    if loc is None:
        return None
    return f"benzofuran-{loc}-amine", f"苯并呋喃-{loc}-胺"


def _btzam_en(pre_en: str, loc: int) -> str:
    stem = f"1,3-benzothiazol-{loc}-amine"
    if not pre_en:
        return stem
    return f"{pre_en}-{stem}" if stem[:1].isdigit() else f"{pre_en}{stem}"


def _btzam_zh(pre_zh: str, loc: int) -> str:
    if not pre_zh:
        return f"{loc}-氨基苯并噻唑"
    return f"{loc}-氨基-{pre_zh}苯并噻唑"


def benzothiazolamine_names(numbered: dict, build_prefix) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant")
    if loc is None:
        return None
    pre_en, pre_zh = build_prefix(numbered.get("substituents") or [], 9, "benzothiazolamine")
    return _btzam_en(pre_en, loc), _btzam_zh(pre_zh, loc)


def _boxam_en(pre_en: str, loc: int) -> str:
    stem = f"1,3-benzoxazol-{loc}-amine"
    if not pre_en:
        return stem
    return f"{pre_en}-{stem}" if stem[:1].isdigit() else f"{pre_en}{stem}"


def _boxam_zh(pre_zh: str, loc: int) -> str:
    if not pre_zh:
        return f"{loc}-氨基苯并噁唑"
    return f"{loc}-氨基-{pre_zh}苯并噁唑"


def benzoxazolamine_names(numbered: dict, build_prefix) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant")
    if loc is None:
        return None
    pre_en, pre_zh = build_prefix(numbered.get("substituents") or [], 9, "benzoxazolamine")
    return _boxam_en(pre_en, loc), _boxam_zh(pre_zh, loc)


def _bimam_en(pre_en: str, loc: int) -> str:
    stem = f"1H-benzimidazol-{loc}-amine"
    if not pre_en:
        return stem
    return f"{pre_en}-{stem}" if stem[:1].isdigit() else f"{pre_en}{stem}"


def _bimam_zh(pre_zh: str, loc: int) -> str:
    stem = f"1H-苯并咪唑-{loc}-胺"
    return f"{pre_zh}-{stem}" if pre_zh else stem


def benzimidazolamine_names(numbered: dict, build_prefix) -> tuple[str, str] | None:
    loc = numbered.get("amine_locant")
    if loc is None:
        return None
    pre_en, pre_zh = build_prefix(numbered.get("substituents") or [], 9, "benzimidazolamine")
    return _bimam_en(pre_en, loc), _bimam_zh(pre_zh, loc)


def benzothiophenol_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("oh_locant")
    if loc is None:
        return None
    return f"1-benzothiophen-{loc}-ol", f"苯并[b]噻吩-{loc}-醇"


_Q_LOCANTS = (1, 2, 3, 4, None, 5, 6, 7, 8, None)


def _q_sub_loc(chain: list[int], attach: int) -> int | None:
    if attach not in chain or len(chain) != 10:
        return None
    loc = _Q_LOCANTS[chain.index(attach)]
    return loc


def _q_cooh_loc(numbered: dict) -> int | None:
    """COOH ring attach with naphthalene-style locants (skip 4a/8a)."""
    parent = numbered.get("parent") or {}
    chain, attach = parent.get("chain") or [], parent.get("ring_attach_idx")
    return None if attach is None else _q_sub_loc(chain, attach)


def quinolinecarboxylic_names(numbered: dict) -> tuple[str, str] | None:
    loc = _q_cooh_loc(numbered)
    if loc is None:
        return None
    return f"quinoline-{loc}-carboxylic acid", f"喹啉-{loc}-甲酸"


def quinolinol_names(numbered: dict) -> tuple[str, str] | None:
    loc = numbered.get("oh_locant")
    if loc is None:
        return None
    return f"quinolin-{loc}-ol", f"喹啉-{loc}-醇"


_IZ_LOCANTS = (1, 2, 3, None, 4, 5, 6, 7, None)


def _iz_fg_loc(numbered: dict) -> int | None:
    """Indazole ring FG attach with indole-style locants (skip 3a/7a)."""
    parent = numbered.get("parent") or {}
    chain, attach = parent.get("chain") or [], parent.get("ring_attach_idx")
    if attach is None or attach not in chain or len(chain) != 9:
        return None
    return _IZ_LOCANTS[chain.index(attach)]


def indazolecarbonitrile_names(numbered: dict) -> tuple[str, str] | None:
    loc = _iz_fg_loc(numbered)
    if loc is None:
        return None
    return f"1H-indazole-{loc}-carbonitrile", f"1H-吲唑-{loc}-甲腈"


def indazolecarbaldehyde_names(numbered: dict) -> tuple[str, str] | None:
    loc = _iz_fg_loc(numbered)
    if loc is None:
        return None
    return f"1H-indazole-{loc}-carbaldehyde", f"1H-吲唑-{loc}-甲醛"


def indolecarboxylic_names(numbered: dict) -> tuple[str, str] | None:
    loc = _iz_fg_loc(numbered)
    if loc is None:
        return None
    return f"1H-indole-{loc}-carboxylic acid", f"1H-吲哚-{loc}-甲酸"


def naphthalenecarboxylic_names(numbered: dict) -> tuple[str, str] | None:
    loc = _q_cooh_loc(numbered)
    if loc is None:
        return None
    return f"naphthalene-{loc}-carboxylic acid", f"萘-{loc}-甲酸"


_PYRIDINE_KIND_FN = {
    "pyridinecarboxylic": pyridinecarboxylic_names,
    "pyridinecarbonitrile": pyridinecarbonitrile_names,
    "pyridinamine": pyridinamine_names,
    "pyridinol": pyridinol_names,
    "benzofuranamine": benzofuranamine_names,
    "benzothiophenol": benzothiophenol_names,
    "quinolinol": quinolinol_names,
    "quinolinecarboxylic": quinolinecarboxylic_names,
    "indazolecarbonitrile": indazolecarbonitrile_names,
    "indazolecarbaldehyde": indazolecarbaldehyde_names,
    "indolecarboxylic": indolecarboxylic_names,
    "naphthalenecarboxylic": naphthalenecarboxylic_names,
}


def pyridine_kind_names(kind: str, numbered: dict, build_prefix) -> tuple[str, str] | None:
    if kind == "pyrimidinamine":
        return pyrimidinamine_names(numbered, build_prefix)
    if kind == "benzothiazolamine":
        return benzothiazolamine_names(numbered, build_prefix)
    if kind == "benzoxazolamine":
        return benzoxazolamine_names(numbered, build_prefix)
    if kind == "benzimidazolamine":
        return benzimidazolamine_names(numbered, build_prefix)
    fn = _PYRIDINE_KIND_FN.get(kind)
    return fn(numbered) if fn is not None else None

