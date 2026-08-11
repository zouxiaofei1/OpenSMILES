"""Retained / multi-substituted benzene parent and prefix helpers (P-22.1.3)."""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
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


def phenyl_parent_names(numbered: dict) -> tuple[str, str]:
    """Free-radical benzene (P-41): pure phenyl stem, prefix comes from _prefix_for.

    No retained names (toluene/anisole/xylene/iso): a phenyl radical is always
    {leaf-locants}{leaves}phenyl, e.g. 4-chlorophenyl, 3,4-dichlorophenyl.
    """
    return "phenyl", "苯基"


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


def _ester_alkoxy_from(numbered) -> tuple[str, str]:
    """O-side alkyl: linear 走 parent.alkoxy_n 保留名表; 特殊基团/复杂回落 o_side 取代基.

    linear (如 hexadecan-16-yl 的 o_side 命名带位次) 必须走 ester_alkyl_en/zh 的
    "hexadecyl"/"十六", 否则长链 O 侧会带 -16-yl 位次 (ester/benzoate 共用此路径).
    """
    if not numbered:
        return "", ""
    parent = numbered.get("parent") or {}
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    if len(o) <= 1 and parent.get("alkoxy_n") is not None:
        from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh
        en, zh = ester_alkyl_en(parent["alkoxy_n"]), ester_alkyl_zh(parent["alkoxy_n"])
        if en and zh:
            return en, zh
    if not o:
        return "", ""
    en0, zh0 = o[0].get("en") or "", (o[0].get("zh") or "").rstrip("基")
    if len(o) == 1:
        return en0, zh0
    # diester: identical arms aggregate (dimethyl / 二甲)
    if all(s.get("en") == en0 for s in o):
        me, mz = MULT_EN.get(len(o), ""), MULT_ZH.get(len(o), "")
        return f"{me}{en0}", f"{mz}{zh0}"
    return en0, zh0


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    """Acid stem + O-side alkyl: methyl butanoate / 丁酸甲酯; stereo after alkyl."""
    en, zh = names
    alk_en, alk_zh = _ester_alkoxy_from(numbered)
    st, body = _stereo_lead(en)
    mid = f"{pre_en}{body}" if pre_en else body
    en = f"{alk_en} {st}{mid}" if alk_en else f"{st}{mid}"
    stz, bodyz = _stereo_lead(zh)
    midz = f"{stz}{pre_zh}{bodyz}" if pre_zh else f"{stz}{bodyz}"
    zh = f"{midz}{alk_zh}酯" if alk_zh else midz
    return en, zh


def join_benzoate_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    """O-side alkyl + benzoate acid stem: ethyl 4-chlorobenzoate / 4-氯苯甲酸乙酯.

    与 ester 同构: 母体只含酸部分, 烷氧基从 o_side 取代基取; zh 恒拼"酯"
    (无烷氧基时保留"苯甲酸酯", 如复杂 O-烷基场景).
    """
    en, zh = names
    alk_en, alk_zh = _ester_alkoxy_from(numbered)
    body = f"{pre_en}{en}" if pre_en else en
    en = f"{alk_en} {body}" if alk_en else body
    midz = f"{pre_zh}{zh}" if pre_zh else zh
    return en, f"{midz}{alk_zh}酯"


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str]:
    if kind in ("ester", "diester"):
        return join_ester_name(pre[0], pre[1], names, numbered)
    if kind == "benzoate":
        return join_benzoate_name(pre[0], pre[1], names, numbered)
    en = join_parent_name(pre[0], names[0])
    zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    return en, zh


def zh_1h_parent(en_parent: str, zh_parent: str, prefix: str) -> str:
    """Prefix Chinese retained 1H-parents with 1H- when ring is substituted."""
    if not prefix or not en_parent.startswith("1H-") or zh_parent.startswith("1H-"):
        return zh_parent
    return f"1H-{zh_parent}"


def benzenediol_names(locs: list[int] | None) -> tuple[str, str] | None:
    if not locs or len(locs) != 2:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"benzene-{loc}-diol", f"苯-{loc}-二酚"




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


# --- Generalized arene FG parent names (P-63.1.4 / P-62.2.1) ---

_ARENE_FG_STEM: dict[str, tuple[str, str, str, str]] = {
    "naphthalenol": ("naphthalen", "萘", "ol", "酚"),
    "naphthalenediol": ("naphthalene", "萘", "diol", "二酚"),
    "pyrazolamine": ("pyrazol", "吡唑", "amine", "胺"),
    "thiazolamine": ("thiazol", "噻唑", "amine", "胺"),
    "quinazolinamine": ("quinazolin", "喹唑啉", "amine", "胺"),
}


def _arene_fg_mono_names(stem_en, stem_zh, suf_en, suf_zh, numbered, kind):
    """Mono-FG arene naming: {stem}-{locant}-{suffix}."""
    loc = numbered.get("oh_locant") if suf_en == "ol" else numbered.get("amine_locant")
    return (f"{stem_en}-{loc}-{suf_en}", f"{stem_zh}-{loc}-{suf_zh}") if loc else None


def _arene_fg_poly_names(stem_en, stem_zh, suf_en, suf_zh, numbered):
    """Poly-FG arene naming: {stem}-{locants}-{suffix}."""
    locs = numbered.get("oh_locants")
    if locs is None: return None
    loc_str = ",".join(str(x) for x in locs)
    return f"{stem_en}-{loc_str}-{suf_en}", f"{stem_zh}-{loc_str}-{suf_zh}"


def _arene_fg_parent_names(kind: str, numbered: dict) -> tuple[str, str] | None:
    """Generalized arene FG parent: {stem}-{locants}-{suffix}."""
    stem = _ARENE_FG_STEM.get(kind)
    if stem is None: return None
    stem_en, stem_zh, suf_en, suf_zh = stem
    if suf_en in ("ol", "amine"):
        return _arene_fg_mono_names(stem_en, stem_zh, suf_en, suf_zh, numbered, kind)
    return _arene_fg_poly_names(stem_en, stem_zh, suf_en, suf_zh, numbered)


_PYRIDINE_KIND_FN = {
    "benzofuranamine": benzofuranamine_names,
    "benzothiophenol": benzothiophenol_names,
    # Generalized arene FG parents
    "naphthalenol": lambda n: _arene_fg_parent_names("naphthalenol", n),
    "naphthalenediol": lambda n: _arene_fg_parent_names("naphthalenediol", n),
    "pyrazolamine": lambda n: _arene_fg_parent_names("pyrazolamine", n),
    "thiazolamine": lambda n: _arene_fg_parent_names("thiazolamine", n),
    "quinazolinamine": lambda n: _arene_fg_parent_names("quinazolinamine", n),
}


def pyridine_kind_names(kind: str, numbered: dict, build_prefix) -> tuple[str, str] | None:
    if kind == "benzothiazolamine":
        return benzothiazolamine_names(numbered, build_prefix)
    if kind == "benzoxazolamine":
        return benzoxazolamine_names(numbered, build_prefix)
    if kind == "benzimidazolamine":
        return benzimidazolamine_names(numbered, build_prefix)
    fn = _PYRIDINE_KIND_FN.get(kind)
    return fn(numbered) if fn is not None else None

