"""L5 fused isocyanato/isothiocyanatobenzene parents (P-61.9 arene)."""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key

_ISO_KINDS = frozenset({"isocyanato", "isothiocyanato"})
_ISO_PARENT = {
    "isocyanato": ("isocyanatobenzene", "异氰酸苯酯"),
    "isothiocyanato": ("isothiocyanatobenzene", "异硫氰酸苯酯"),
}
_ISO_FUSED_OK = frozenset({"alkyl", "trifluoromethyl"})


def _iso_subs(subs: list) -> list:
    return [s for s in subs if s.get("kind") in _ISO_KINDS]


def is_fused_iso(numbered: dict) -> bool:
    """Mono iso + only alkyl/CF3 → fuse iso into parent (P-61.9 arene)."""
    subs = numbered.get("substituents") or []
    if len(_iso_subs(subs)) != 1:
        return False
    return all(s.get("kind") in _ISO_KINDS | _ISO_FUSED_OK for s in subs)


def _iso_reloc(loc: int, iso: int, rev: bool) -> int:
    d = (iso - loc) % 6 if rev else (loc - iso) % 6
    return d + 1


def _iso_dir_key(others: list, iso: int, rev: bool) -> tuple:
    ordered = sorted(
        (s for s in others if "locant" in s),
        key=lambda s: alkyl_alpha_key(s.get("en") or ""),
    )
    return tuple(_iso_reloc(int(s["locant"]), iso, rev) for s in ordered)


def _iso_pick_rev(others: list, iso: int) -> bool:
    return _iso_dir_key(others, iso, True) < _iso_dir_key(others, iso, False)


def _iso_renum(others: list, iso_loc: int, rev: bool) -> list:
    return [{**s, "locant": _iso_reloc(int(s["locant"]), iso_loc, rev)} for s in others]


def _iso_others_renum(subs: list, iso: dict) -> list:
    loc = int(iso["locant"])
    others = [s for s in subs if s is not iso]
    return _iso_renum(others, loc, _iso_pick_rev(others, loc))


def _iso_en_paren_cf3(en: str) -> str:
    """Gold EN uses 4-(trifluoromethyl)...; ZH omits parens."""
    if "trifluoromethyl" not in en or "(trifluoromethyl)" in en:
        return en
    return en.replace("trifluoromethyl", "(trifluoromethyl)")


def iso_fused_parent(numbered: dict) -> tuple[str, str]:
    kind = _iso_subs(numbered.get("substituents") or [])[0].get("kind")
    return _ISO_PARENT[kind]


def iso_fused_prefix(numbered: dict, build_prefix) -> tuple[str, str]:
    """Keep relative locants; iso is 1 (no mono-benzene omit)."""
    subs = numbered.get("substituents") or []
    iso = _iso_subs(subs)[0]
    renum = _iso_others_renum(subs, iso)
    if not renum:
        return "", ""
    en, zh = build_prefix(renum, 6, None)
    return _iso_en_paren_cf3(en), zh
