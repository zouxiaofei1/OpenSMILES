"""L5 取代基前缀分组与双语渲染。"""
from __future__ import annotations

from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.constants import MULT_EN, MULT_ZH

def _group_by_stem(substituents: list) -> dict[str, list]:
    """按 en 词干对取代基分组，返回词干到列表的映射。"""
    groups: dict[str, list] = {}
    for s in substituents:
        groups.setdefault(s.get("en") or "", []).append(s)
    return groups


def _locant_str(subs: list) -> str:
    """对位次排序并拼接成逗号分隔串（如 1,3）。"""
    locs = sorted(int(s["locant"]) for s in subs if "locant" in s)
    return ",".join(str(x) for x in locs)


def _mult_en(n: int) -> str:
    """英文数量前缀（di/tri…），查表返回。"""
    return MULT_EN.get(n, "")


def _mult_zh(n: int) -> str:
    """中文数量前缀（二/三…），查表返回。"""
    return MULT_ZH.get(n, "")


_KEEP_LOCANT_KINDS = frozenset({
    "acid",
    "acetophenone", "benzoyl_chloride", "benzoyl_bromide",
})


def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None,
                      scaffold: str | None = None, has_ene: bool = False) -> bool:
    """判断取代基位次可否省略（环烷烃/苯单取代、酰胺 N- 等情形）。"""
    if kind == "phenyl":
        # 苯基自由基：连接点隐含为 locant 1，因此每个叶子保留其位次（4-chlorophenyl，而非 chlorophenyl）。
        return False
    # 纯烃环单取代位次隐含：环烷烃/苯 base 的 kind 均收敛为 alkane，环系由 scaffold_id 承载；环烯取代基位次必须保留（1-methylcyclohexene）。
    if n_carbons <= 1 or (
        kind == "alkane" and scaffold in ("carbocycle", "benzene") and not has_ene
    ) and len(substituents) == 1:
        return True
    if kind == "amide":
        return {s.get("kind") for s in substituents} <= {
            "n_alkyl", "n_phenyl", "n_benzyl", "n_block",
        }
    if kind in _KEEP_LOCANT_KINDS or kind == "ketone":
        return False
    if any(s.get("paren") or (s.get("en") or "")[:1] == "(" for s in substituents):
        return False
    return n_carbons == 2 and len(substituents) == 1


def _stem_needs_paren(stem: str, subs: list, omit: bool) -> bool:
    """加括号：显式标记、前导位次词干，或多个 CF3（英文）。"""
    if any(s.get("paren") for s in subs):
        return True
    if stem and stem[0].isdigit():
        return True
    return (not omit) and stem == "trifluoromethyl"


def _wrap_stem(stem: str, need: bool) -> str:
    """按 need 给词干加括号（词干含括号时改用方括号）。"""
    return stem if not need else (f"[{stem}]" if "(" in stem else f"({stem})")


def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    """拼单个英文前缀：数量 + 词干（可省略位次时省略 locant）。"""
    mult = _complex_mult_en(stem, len(subs)) or _mult_en(len(subs))
    s = _wrap_stem(stem, _stem_needs_paren(stem, subs, omit))
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"


def _complex_mult_en(stem: str, n: int) -> str:
    """英文复杂数量前缀：carboxy 词干用 bis/tris。"""
    return {2: "bis", 3: "tris"}.get(n, "") if "carboxy" in stem else ""


def _complex_mult_zh(stem: str, n: int) -> str:
    """中文复杂数量前缀：羧基词干用 双/三。"""
    return {2: "双", 3: "三"}.get(n, "") if "羧" in stem else ""


def _prefix_one_zh(zh_stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> str:
    """拼单个中文前缀：数量 + 词干（含 CF3 括号特例）。"""
    mult = _complex_mult_zh(zh_stem, len(subs)) or _mult_zh(len(subs))
    en = subs[0].get("en") or ""
    need = any(s.get("paren") for s in subs) or (en[:1].isdigit() if en else False)
    if paren_cf3 and zh_stem == "三氟甲基":
        need = True
    s = _wrap_stem(zh_stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"


def _sorted_stems(groups: dict[str, list]) -> list[str]:
    """按烷基字母键排序非空词干列表。"""
    return sorted((k for k in groups if k), key=alkyl_alpha_key)


def _parts_for_stem(stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> tuple[str, str]:
    """按词干生成中英文前缀（N- 类取代基强制省略位次）。"""
    zh_stem = subs[0].get("zh") or ""
    if (subs[0].get("kind") or "") in ("n_alkyl", "n_phenyl", "n_benzyl", "n_block"):
        omit = True
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit, paren_cf3)


def _collect_parts(groups: dict[str, list], omit: bool, paren_cf3: bool = False) -> tuple[list[str], list[str]]:
    """汇总所有词干的中英文前缀部件列表。"""
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for stem in _sorted_stems(groups):
        en_p, zh_p = _parts_for_stem(stem, groups[stem], omit, paren_cf3)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts


def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None,
                  scaffold: str | None = None, has_ene: bool = False) -> tuple[str, str]:
    """组合完整取代基前缀：滤 O 侧、判 omit、按词干分组拼接。"""
    if not substituents:
        return "", ""
    # ester 的 O 侧烷基由 join_kind_name 作为烷氧基臂消费，永不作前缀
    substituents = [s for s in substituents if not s.get("o_side")]
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents, kind, scaffold, has_ene)
    paren = scaffold == "benzene" and kind == "alkane" and len(substituents) >= 4
    en_parts, zh_parts = _collect_parts(_group_by_stem(substituents), omit, paren)
    return "-".join(en_parts), "-".join(zh_parts)

def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    """从 numbered 提取母体上下文并委托 _build_prefix 构建前缀。"""
    parent = numbered.get("parent") or {}
    has_ene = bool(parent.get("double_bond") or parent.get("double_bonds"))
    return _build_prefix(numbered.get("substituents") or [], n, kind,
                         parent.get("scaffold_id"), has_ene)
