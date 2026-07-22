"""Layer 5 substituent-prefix grouping and bilingual rendering."""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.layer5.benzene_names import benzene_prefix
from namepredict.constants import MULT_EN, MULT_ZH
_H5COOH_KINDS = frozenset({"furancarboxylic", "thiophenecarboxylic", "pyrrolecarboxylic", "imidazolecarboxylic", "pyrazolecarboxylic"})
_SHCOOH_KINDS = frozenset({"piperidinecarboxylic", "pyrrolidinecarboxylic", "piperazinecarboxylic", "morpholinecarboxylic", "oxolanecarboxylic", "oxanecarboxylic", "thiolanecarboxylic", "aziridinecarboxylic"})

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
    "benzoic", "benzaldehyde", "acetophenone", "pyridinecarboxylic",
    "pyridinecarbonitrile", "benzoate", "benzonitrile", "benzoyl_chloride",
    "benzoyl_bromide", "benzamide", "cycloalkanecarboxylic", "cycloalkanecarbaldehyde",
    "cycloalkanecarbonitrile", "cycloalkanecarboxamide",
    "cycloalkanecarboxylate", "cycloalkanecarbonyl_chloride",
    "cycloalkanecarbonyl_bromide"}) | _H5COOH_KINDS | _SHCOOH_KINDS
def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None) -> bool:
    if n_carbons <= 1 or (kind in ("cycloalkane", "benzene") and len(substituents) == 1):
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
def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None) -> tuple[str, str]:
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents, kind)
    paren = kind == "benzene" and len(substituents) >= 4
    en_parts, zh_parts = _collect_parts(_group_by_stem(substituents), omit, paren)
    return "-".join(en_parts), "-".join(zh_parts)
def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    if kind == "benzene":
        return benzene_prefix(numbered, _build_prefix)
    _skip = ("benzoate", "pyrimidinamine", "benzothiazolamine",
             "benzoxazolamine", "benzimidazolamine", "boronic")
    if kind in _skip:
        return "", ""
    from namepredict.layer5.acyl_halide_names import _is_isobutyryl
    if kind == "acyl_bromide" and _is_isobutyryl(numbered):
        return "", ""
    return _build_prefix(numbered.get("substituents") or [], n, kind)
