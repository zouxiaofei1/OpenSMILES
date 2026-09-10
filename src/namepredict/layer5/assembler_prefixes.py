"""L5 取代基前缀分组与双语渲染。"""
from __future__ import annotations

import re

from namepredict.layer1 import fg_registry as _fg_reg
from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.constants import MULT_EN, MULT_ZH, N_PREFIX_KINDS

def _group_by_stem(substituents: list) -> dict[str, list]:
    """按 en 词干对取代基分组，返回词干到列表的映射。"""
    groups: dict[str, list] = {}
    for s in substituents:
        groups.setdefault(s.get("en") or "", []).append(s)
    return groups


def _locant_str(subs: list) -> str:
    """按位次排序拼接成逗号串；N-型取代基渲染为字母位次 N（与 C 数字位次并排，N 自然排最前）。"""
    from namepredict.layer4.locant_key import locant_str_sort
    tokens = []
    for s in subs:
        if "locant" not in s:
            continue
        kind = s.get("kind") or ""
        tokens.append("N" if kind in N_PREFIX_KINDS else s["locant"])
    return ",".join(str(x) for x in locant_str_sort(tokens))


def _mult_en(n: int) -> str:
    """英文数量前缀（di/tri…），查表返回。"""
    return MULT_EN.get(n, "")


def _mult_zh(n: int) -> str:
    """中文数量前缀（二/三…），查表返回。"""
    return MULT_ZH.get(n, "")


_KEEP_LOCANT_KINDS = _fg_reg.keep_locant_fgs()


def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None,
                      scaffold: str | None = None, has_ene: bool = False) -> bool:
    """判断取代基位次可否省略（环烷烃/苯单取代、酰胺 N- 等情形）。"""
    if kind == "radical":  # 自由基母体：连接点隐含为 locant 1，每个叶子保留位次（4-chlorophenyl）；单碳链（methyl 型）取代基都在唯一 C1，位次无信息量故省略（(phenyloxy)methyl）。
        return n_carbons == 1
    if n_carbons <= 1:  # 单碳母体位次隐含省略；但 N- 型与 C- 型取代基共存时 C 侧必须带数字位次消歧（P-62.2.4.1.2：胺的数字位次含单核母体的 '1'，与 'N' 位次并引，1,1-dimethoxy-N,N-dimethylmethanamine）。
        kinds = {(s.get("kind") or "") for s in substituents}
        return not (kinds & N_PREFIX_KINDS and kinds - N_PREFIX_KINDS)
    if (  # 纯烃环单取代位次隐含：环烷烃/苯 base 的 kind 均收敛为 alkane，环系由 scaffold_id 承载；环烯取代基位次必须保留（1-methylcyclohexene）。
        kind == "alkane" and scaffold in ("carbocycle", "benzene") and not has_ene
    ) and len(substituents) == 1:
        return True
    if kind == "amide":
        return {s.get("kind") for s in substituents} <= N_PREFIX_KINDS
    if kind in _KEEP_LOCANT_KINDS:
        return False
    if any(s.get("paren") or (s.get("en") or "")[:1] == "(" for s in substituents):
        return False
    complex_sub = any(  # 复合取代基（自带位次如 1H-indol-5-yl / propan-2-ylsulfanyl，或显式括号）须保留母体 2- 消歧；
        s.get("paren") or (s.get("en") or "")[:1] == "("        # 简单 FG 前缀（amino/hydroxy/chloro）位次无信息量，省略。
        or re.search(r"\d", s.get("en") or "") for s in substituents
    )
    return (  # C2 单取代省略仅对端碳(FG 所在 C1)无可取代 H 的母体成立(腈/酯/醛/酰胺等)；醇/胺/硫醇的 C1 带可取代 H，2- 位取代构成不同异构体(P-14.3.4.4)，2- 必须保留。
        n_carbons == 2 and len(substituents) == 1
        and kind not in ("alcohol", "amine", "thiol")
        and not complex_sub
    )


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


# O/S 桥后缀（gold 平铺式 -yl]oxy/-yl]sulfanyl：括号闭在 -yl 后、后缀放括号外，见 P-63.2.2）。
_BRIDGE_SUFFIX_EN = ("oxy", "sulfanyl")


def _split_bridge_suffix(stem: str) -> tuple[str, str] | None:
    """带立体描述符的基 -<N>-yl]oxy/-yl]sulfanyl 拆分：(base-yl, 桥后缀)；仅拆 base 含手性描述符(如 2R/3S)的情形，无手性 acyclic/苄基(…methylsulfanyl/…propan-2-yloxy 等)整括不拆。"""
    for suf in _BRIDGE_SUFFIX_EN:
        if not stem.endswith(suf):
            continue
        base = stem[: -len(suf)]
        if not (base.endswith("yl") and re.search(r"\d[RrSs]", base)):
            continue
        return base, suf
    return None


def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    """拼单个英文前缀：数量 + 词干（可省略位次时省略 locant）。"""
    mult = _complex_mult_en(stem, subs, len(subs)) or _mult_en(len(subs))
    need = _stem_needs_paren(stem, subs, omit)
    if not omit and not mult and need:
        sp = _split_bridge_suffix(stem)
        if sp is not None:  # O/S 桥平铺式：括号闭在 -yl 后，-oxy/-sulfanyl 追加在括号外（gold 449:0 形式）。
            base, suf = sp
            return f"{_locant_str(subs)}-{_wrap_stem(base, True)}{suf}"
    s = _wrap_stem(stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"


def _is_compound_mult(stem: str, subs: list) -> bool:
    """待倍增组分是否为复合/被取代前缀(须用 bis/tris 而非 di/tri, P-16.3.2)：retained 组合叶(carboxy-/hydroxymethyl 等整叶名含修饰前缀)由词干子串兜底，递归命名/括号组分由组成员 paren 标记体现。"""
    return ("carboxy" in stem) or any(s.get("paren") for s in subs)


def _complex_mult_en(stem: str, subs: list, n: int) -> str:
    """英文复杂数量前缀：复合组分用 bis/tris/tetrakis。"""
    return {2: "bis", 3: "tris", 4: "tetrakis"}.get(n, "") if _is_compound_mult(stem, subs) else ""


def _complex_mult_zh(stem: str, subs: list, n: int) -> str:
    """中文复杂数量前缀：复合组分用 双/三/四。"""
    return {2: "双", 3: "三", 4: "四"}.get(n, "") if (("羧" in stem) or any(s.get("paren") for s in subs)) else ""


def _prefix_one_zh(zh_stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> str:
    """拼单个中文前缀：数量 + 词干（CF3 特例：简单氟代甲基不加括号）。"""
    mult = _complex_mult_zh(zh_stem, subs, len(subs)) or _mult_zh(len(subs))
    en = subs[0].get("en") or ""
    need = any(s.get("paren") for s in subs) or (en[:1].isdigit() if en else False)  # 停用：zh_stem == "三氟甲基" 时 need = False（简单氟代甲基不加括号）
    s = _wrap_stem(zh_stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"


def _sorted_stems(groups: dict[str, list]) -> list[str]:
    """按烷基字母键排序非空词干列表。"""
    return sorted((k for k in groups if k), key=alkyl_alpha_key)


def _n_prefix_en(n: int, stem: str) -> str:
    """英文 N- 前缀：N-methyl / N,N-dimethyl / N,N,N-trimethyl。"""
    if n == 1:
        return f"N-{stem}"
    ns = ",".join("N" for _ in range(n))
    return f"{ns}-{_mult_en(n)}{stem}"


def _n_prefix_zh(n: int, stem: str) -> str:
    """中文 N- 前缀：N-甲基 / N,N-二甲基。"""
    if n == 1:
        return f"N-{stem}"
    ns = ",".join("N" for _ in range(n))
    return f"{ns}-{_mult_zh(n)}{stem}"


def _parts_for_stem(stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> tuple[str, str]:
    """按词干生成中英文前缀（N- 类取代基加 N- 前缀并强制省略位次）。"""
    zh_stem = subs[0].get("zh") or ""
    if subs and all((s.get("kind") or "") in N_PREFIX_KINDS for s in subs):  # 整组全为 N-型取代基才走 N-计数前缀；同词干混入 C-型时落入数字通道，N-型成员由 _locant_str 渲染为 N（如 N,N,2-trimethyl）。
        need = any(s.get("paren") for s in subs) or bool(stem and stem[0].isdigit())  # 复合取代基（含 locant 位次/显式 paren）须整体加括号：N-(3-bromophenyl)。
        s_en = _wrap_stem(stem, need)
        s_zh = _wrap_stem(zh_stem, need)
        return _n_prefix_en(len(subs), s_en), _n_prefix_zh(len(subs), s_zh)
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit, paren_cf3)


def _collect_parts(groups: dict[str, list], omit: bool, paren_cf3: bool = False,
                   bracket: bool = False) -> tuple[list[str], list[str]]:
    """汇总所有词干的中英文前缀部件列表；bracket(单碳多不同取代, P-16.5.1.3.1)时第二词干起整体加圆括号(倍增前缀不括入)，词干间无连字符。"""
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for i, stem in enumerate(_sorted_stems(groups)):
        subs = groups[stem]
        if bracket and i >= 1:
            en_parts.append(f"{_mult_en(len(subs))}({stem})")
            zh_parts.append(f"{_mult_zh(len(subs))}({subs[0].get('zh') or ''})")
            continue
        en_p, zh_p = _parts_for_stem(stem, subs, omit, paren_cf3)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts


def _groups_simple(groups: dict[str, list]) -> bool:
    """全部词干为简单取代基(无显式括号、无数字/locant 前导、非 N- 类)才适用括号式。"""
    for subs in groups.values():
        for s in subs:
            if s.get("paren") or (s.get("en") or "")[:1].isdigit():
                return False
            if (s.get("kind") or "") in N_PREFIX_KINDS:
                return False
    return True


def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None,
                  scaffold: str | None = None, has_ene: bool = False) -> tuple[str, str]:
    """组合完整取代基前缀：滤 O 侧、判 omit、按词干分组拼接。"""
    if not substituents:
        return "", ""
    substituents = [s for s in substituents if not s.get("o_side")]  # ester 的 O 侧烷基由 join_kind_name 作为烷氧基臂消费，永不作前缀
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents, kind, scaffold, has_ene)
    paren = scaffold == "benzene" and kind == "alkane" and len(substituents) >= 4
    groups = _group_by_stem(substituents)
    bracket = bool(omit) and n_carbons == 1 and kind == "radical" \
        and len(groups) >= 2 and _groups_simple(groups)
    en_parts, zh_parts = _collect_parts(groups, omit, paren, bracket)  # P-16.5.1.3.1/.3.2：单碳(meth)母链带 ≥2 个不同简单取代基且位次省略 → 首词干平铺、第二及以后各自括号；单碳链取代基必同处唯一碳，括号式即 locant 省略时的消歧写法。
    sep = "" if bracket else "-"
    return sep.join(en_parts), sep.join(zh_parts)

def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    """从 numbered 提取母体上下文并委托 _build_prefix 构建前缀。"""
    parent = numbered.get("parent") or {}
    has_ene = bool(parent.get("double_bond") or parent.get("double_bonds"))
    if kind == "radical" and parent.get("radical_anchor_element"):  # 杂原子锚点自由基：烷基取代基已并入组装名（ethyloxy），不再加前缀。
        return "", ""
    return _build_prefix(numbered.get("substituents") or [], n, kind,
                         parent.get("scaffold_id"), has_ene)
