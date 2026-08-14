"""Layer 5 substituent-prefix grouping and bilingual rendering."""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.layer5.benzene_names import benzene_prefix
from namepredict.constants import MULT_EN, MULT_ZH
def _group_by_stem(substituents: list) -> dict[str, list]:
    groups: dict[str, list] = {}
    for s in substituents:
        groups.setdefault(s.get("en") or "", []).append(s)
    return groups
def _locant_str(subs: list) -> str:
    locs = sorted(int(s["locant"]) for s in subs if "locant" in s)
    return ",".join(str(x) for x in locs)
def _mult_en(n: int) -> str: return MULT_EN.get(n, "")
def _mult_zh(n: int) -> str: return MULT_ZH.get(n, "")
_KEEP_LOCANT_KINDS = frozenset({
    "acid",
    "benzoic", "benzaldehyde", "acetophenone", "benzoate", "benzonitrile",
    "benzoyl_chloride", "benzoyl_bromide", "benzamide",
})
def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None,
                      scaffold: str | None = None, has_ene: bool = False) -> bool:
    if kind == "phenyl":
        # Phenyl radical: attach is implicit locant 1, so every leaf keeps its
        # locant (4-chlorophenyl, not chlorophenyl).
        return False
    # 纯烃饱和环单取代（cycloalkane，kind 收敛为 alkane + carbocycle scaffold）位次隐含；
    # 环烯取代基位次必须保留（1-methylcyclohexene）。
    if n_carbons <= 1 or (
        (kind == "alkane" and scaffold == "carbocycle" and not has_ene) or kind == "benzene"
    ) and len(substituents) == 1:
        return True
    if kind in ("sec_amine", "tert_amine", "amide", "benzamide"):
        return {s.get("kind") for s in substituents} <= {
            "n_alkyl", "n_phenyl", "n_benzyl", "n_block",
        }
    if kind in _KEEP_LOCANT_KINDS or kind == "ketone":
        return False
    if any(s.get("paren") or (s.get("en") or "")[:1] == "(" for s in substituents):
        return False
    return n_carbons == 2 and len(substituents) == 1
def _stem_needs_paren(stem: str, subs: list, omit: bool) -> bool:
    """Paren: explicit flag, leading-locant stems, or multi CF3 (EN)."""
    if any(s.get("paren") for s in subs):
        return True
    if stem and stem[0].isdigit():
        return True
    return (not omit) and stem == "trifluoromethyl"
def _wrap_stem(stem: str, need: bool) -> str: return stem if not need else (f"[{stem}]" if "(" in stem else f"({stem})")
def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    mult = _complex_mult_en(stem, len(subs)) or _mult_en(len(subs))
    s = _wrap_stem(stem, _stem_needs_paren(stem, subs, omit))
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"
def _complex_mult_en(stem: str, n: int) -> str:
    return {2: "bis", 3: "tris"}.get(n, "") if "carboxy" in stem else ""
def _complex_mult_zh(stem: str, n: int) -> str:
    return {2: "双", 3: "三"}.get(n, "") if "羧" in stem else ""
def _prefix_one_zh(zh_stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> str:
    mult = _complex_mult_zh(zh_stem, len(subs)) or _mult_zh(len(subs))
    en = subs[0].get("en") or ""
    need = any(s.get("paren") for s in subs) or (en[:1].isdigit() if en else False)
    if paren_cf3 and zh_stem == "三氟甲基":
        need = True
    s = _wrap_stem(zh_stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"
def _sorted_stems(groups: dict[str, list]) -> list[str]: return sorted((k for k in groups if k), key=alkyl_alpha_key)
def _parts_for_stem(stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> tuple[str, str]:
    zh_stem = subs[0].get("zh") or ""
    if (subs[0].get("kind") or "") in ("n_alkyl", "n_phenyl", "n_benzyl", "n_block"):
        omit = True
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit, paren_cf3)
def _collect_parts(groups: dict[str, list], omit: bool, paren_cf3: bool = False) -> tuple[list[str], list[str]]:
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for stem in _sorted_stems(groups):
        en_p, zh_p = _parts_for_stem(stem, groups[stem], omit, paren_cf3)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts
def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None,
                  scaffold: str | None = None, has_ene: bool = False) -> tuple[str, str]:
    if not substituents:
        return "", ""
    # ester O-side alkyl is consumed by join_kind_name as the alkoxy arm, never a prefix
    substituents = [s for s in substituents if not s.get("o_side")]
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents, kind, scaffold, has_ene)
    paren = kind == "benzene" and len(substituents) >= 4
    en_parts, zh_parts = _collect_parts(_group_by_stem(substituents), omit, paren)
    return "-".join(en_parts), "-".join(zh_parts)

def _is_isobutyryl(numbered: dict) -> bool:
    """C3 acyl + single methyl at locant 2 (isobutyryl retained topology)."""
    if int((numbered.get("parent") or {}).get("n_carbons") or 0) != 3:
        return False
    subs = numbered.get("substituents") or []
    if len(subs) != 1:
        return False
    s = subs[0]
    return s.get("en") == "methyl" and int(s.get("locant") or 0) == 2

def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    if kind == "benzene":
        return benzene_prefix(numbered, _build_prefix)
    if kind == "acyl_bromide" and _is_isobutyryl(numbered):
        return "", ""
    parent = numbered.get("parent") or {}
    has_ene = bool(parent.get("double_bond") or parent.get("double_bonds"))
    return _build_prefix(numbered.get("substituents") or [], n, kind,
                         parent.get("scaffold_id"), has_ene)
