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
    """加括号：显式标记、前导位次词干、前导立体描述符，或多个 CF3（英文）。"""
    if any(s.get("paren") for s in subs):
        return True
    if stem and stem[0].isdigit():
        return True
    return (not omit) and (stem == "trifluoromethyl" or bool(_STEREO_LEAD_RE.match(stem)))


def _wrap_stem(stem: str, need: bool) -> str:
    """按 need 给词干加括号（词干含括号时改用方括号）。"""
    return stem if not need else (f"[{stem}]" if "(" in stem else f"({stem})")

_STEREO_LEAD_RE = re.compile(r"\(\d+[RSEZ](?:,\d+[RSEZ])*\)-")  # 取代基名以立体描述符开头：(1Z)-、(2R,4R)-、(9Z,12Z)-。

_BRIDGE_SUFFIX_EN = ("oxy", "sulfanyl", "amino")  # O/S/N 桥后缀（gold 平铺式 -yl]oxy/-yl]amino：括号闭在 -yl 后、后缀放括号外，见 P-63.2.2.1）。

_CHAIN_STEM = (r"(?:meth|eth|prop|but|pent|hex|hept|oct|non|dec|undec|dodec|tridec|tetradec|"
               r"pentadec|hexadec|heptadec|octadec|nonadec|eicos)")
_SIMPLE_CHAIN_YL_RE = re.compile(r"^(?:\d+-)?" + _CHAIN_STEM + r"a?n-\d+-yl$")  # 无取代直链 -yl（propan-2-yl）：gold 与 O/S 桥融合平铺（propan-2-yloxy），不拆围栏。
_SUBST_CHAIN_YL_RE = re.compile(_CHAIN_STEM + r"a?n-\d+-yl$")  # 带取代基的直链 -yl（1,3-dihydroxypropan-2-yl）：gold 仍与桥融合平铺。
_TERMINAL_CHAIN_YL_RE = re.compile(_CHAIN_STEM + r"yl$")  # 自由价在端碳的直链基（…phenyl)methyl、5-(X)pentyl）：前端括号只属于其取代基，与桥融合平铺。
_BENZYL_TAIL_RE = re.compile(r"\]methyl$")  # 苄基型前端（…oxolan-2-yl]methyl）：桥后缀直接缀在甲基上（methylsulfanyl），不拆。


def _front_needs_enclosure(base: str, suf: str) -> bool:
    """O/S/N 桥前端是否为自带围栏的复合取代基：带立体描述符的自由价碳，或（oxy/sulfanyl 桥时）环型内嵌位次 …oxan-2-yl。"""
    if base.endswith(("sulfonyl", "sulfinyl")):  # 磺酰基/亚磺酰基前端的围栏由 L3 定形（复合前端括起、简单前端平铺），此处不再拆。
        return False
    if _BENZYL_TAIL_RE.search(base):  # 苄基型前端（…yl]methylsulfanyl）：桥后缀直接缀在甲基上，不拆。
        return False
    if  re.match(r"^\(\d+[RrSs]", base) and not base.endswith("oyl"):  # 自由价碳带手性描述符（(2R)-2-amino-2-carboxyethyl）：前端必须括起；酰基前端按 …oyloxy 融合（P-63.2.2.1.1 的 acetyloxy/benzoyloxy）。
        return True
    if "[" in base:  # 前端已用方括号。amino 桥细分：括号后仍接数字位次前缀（…]sulfanylethylamino]-3-oxopropyl）时拆——前端是被位次取代的支链，须整体括起（gold 取 [[X]amino]）；括号后接桥/链续写（…enoyl]sulfanylethyl）时 gold 平铺为 sulfanylethylamino，不拆。oxy/sulfanyl 桥不细分（…]phenylsulfanyl、…]acetyloxy 的围栏由 L3 定形）。
        return True if suf != "amino" else bool(
            re.search(r"\]-?\d", base)                            # 括号后接数字位次前缀：…]sulfanylethylamino]-3-oxopropyl
            or re.search(r"-\d+-\[", base[: base.find("[") + 1]))  # 括号前已有数字位次前缀：3-oxo-3-[X]propyl
    if "(" in base:  # 前端自带括号（取代基/立体描述符）；自由价在端碳的链基（benzyl/5-(X)pentyl）平铺。
        return not _TERMINAL_CHAIN_YL_RE.search(base)
    if suf in ("oxy", "sulfanyl") and re.search(r"\d", base) and base.endswith("phenyl"):  #
        return True
    if suf == "amino":  # P-63.2.2.1.2：amino 桥按 HS- 取代式（naphthalen-2-ylamino）融合，无括号前端不拆。
        return False
    return bool(re.search(r"-\d+-yl$", base)) and not _SUBST_CHAIN_YL_RE.search(base)


def _split_bridge_suffix(stem: str) -> tuple[str, str] | None:
    """拆 -yl]o xy/-yl]sulfanyl/-yl]amino 平铺式：(前端, 桥后缀)；括号闭在前端 -yl 后、桥后缀留在括号外。前端为简单保留基（methyl/benzyl）、直链 -yl（propan-2-yl）或酰基时整括不拆（P-63.2.2.1.1）。"""
    for suf in _BRIDGE_SUFFIX_EN:
        if not stem.endswith(suf):
            continue
        base = stem[: -len(suf)]
        if not base.endswith("yl") or not _front_needs_enclosure(base, suf):
            continue
        return base, suf
    return None


def _sbridge_flat_stem(stem: str) -> bool:
    """磺酰基/亚磺酰基桥 + 直链 -yl 前端：英文侧平铺不加围栏（propan-2-ylsulfonyl，P-63.2.1；中文侧仍括注）。"""
    return any(stem.endswith(suf) and _SIMPLE_CHAIN_YL_RE.match(stem[: -len(suf)])
               for suf in ("sulfonyl", "sulfinyl"))


def _bridge_body(base: str, suf: str) -> str:
    """O/S/N 桥平铺式主体：前端自带围栏 + 桥后缀留括号外；前端围栏已是方括号且桥为氨基时整体再括一层（P-63.2.2.1.2：[[X]amino]propanoyl）。"""
    body = f"{_wrap_stem(base, True)}{suf}"
    return f"[{body}]" if body.startswith("[") and suf in ("amino", "氨基") else body


def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    """拼单个英文前缀：数量 + 词干（可省略位次时省略 locant）。"""
    mult = _complex_mult_en(stem, subs, len(subs)) or _mult_en(len(subs))
    need = _stem_needs_paren(stem, subs, omit) and not _sbridge_flat_stem(stem)
    if need:
        sp = _split_bridge_suffix(stem)
        if sp is not None:  # O/S/N 桥平铺式：括号闭在前端 -yl 后，-oxy/-sulfanyl/-amino 追加在括号外（P-63.2.2.1.1）。
            base, suf = sp
            body = _bridge_body(base, suf)
            return f"{mult}{body}" if omit else f"{_locant_str(subs)}-{mult}{body}"
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


def _split_bridge_suffix_zh(zh_stem: str, en_stem: str) -> tuple[str, str] | None:
    """中文侧 O/S/N 桥平铺式拆分（…基]氧基/硫基/氨基）：判据与英文侧同步，仅当英文 stem 拆时才拆，保证中英围栏同形。"""
    if _split_bridge_suffix(en_stem) is None:
        return None
    for suf in ("氧基", "硫基", "氨基"):
        if zh_stem.endswith(suf):
            base = zh_stem[: -len(suf)]
            if base.endswith("基"):
                return base, suf
            if base.endswith("-"):  # 环/链自由价位次在桥融合时被氧基顶掉「基」（喹啉-8-氧基 → 喹啉-8-基），拆时补回。
                return f"{base}基", suf
    return None


def _prefix_one_zh(zh_stem: str, subs: list, omit: bool, paren_cf3: bool = False,
                   en_stem: str = "") -> str:
    """拼单个中文前缀：数量 + 词干（CF3 特例：简单氟代甲基不加括号）。"""
    mult = _complex_mult_zh(zh_stem, subs, len(subs)) or _mult_zh(len(subs))
    en = subs[0].get("en") or ""
    need = any(s.get("paren") for s in subs) or (en[:1].isdigit() if en else False)  # 停用：zh_stem == "三氟甲基" 时 need = False（简单氟代甲基不加括号）
    if not omit and _STEREO_LEAD_RE.match(en):  # 与英文侧同步：前导立体描述符 + 位次须整体围栏（5-[(1Z)-丙-1-烯基]苯）
        need = True
    if need:
        sp = _split_bridge_suffix_zh(zh_stem, en_stem)
        if sp is not None:  # 与英文侧同形：括号闭在前端「基」后，氧基/硫基/氨基留括号外（P-63.2.2.1.1）
            base, suf = sp
            body = _bridge_body(base, suf)
            return f"{mult}{body}" if omit else f"{_locant_str(subs)}-{mult}{body}"
    s = _wrap_stem(zh_stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"


def _sorted_stems(groups: dict[str, list]) -> list[str]:
    """按烷基字母键排序非空词干列表。"""
    return sorted((k for k in groups if k), key=alkyl_alpha_key)


def _n_prime_tokens(subs: list, primes: dict[int, int] | None) -> list[str]:
    """N-型取代基的位次记号列表：同一 N 上为 N，不同 N 上依次 N、N'（P-14.5 多氮位次消歧）。"""
    return sorted(("N" + "'" * ((primes or {}).get(s.get("attach_idx"), 0))) for s in subs)


def _n_prefix_en(n: int, stem: str, tokens: list[str], subs: list) -> str:
    """英文 N- 前缀：N-methyl / N,N-dimethyl / N,N'-bis[...]。"""
    if n == 1:
        return f"{tokens[0]}-{stem}"
    ns = ",".join(tokens)
    return f"{ns}-{_complex_mult_en(stem, subs, n) or _mult_en(n)}{stem}"


def _n_prefix_zh(n: int, zh_stem: str, tokens: list[str], subs: list) -> str:
    """中文 N- 前缀：N-甲基 / N,N-二甲基 / N,N'-双[...]。"""
    if n == 1:
        return f"{tokens[0]}-{zh_stem}"
    ns = ",".join(tokens)
    return f"{ns}-{_complex_mult_zh(zh_stem, subs, n) or _mult_zh(n)}{zh_stem}"


def _parts_for_stem(stem: str, subs: list, omit: bool, paren_cf3: bool = False,
                    primes: dict[int, int] | None = None) -> tuple[str, str]:
    """按词干生成中英文前缀（N- 类取代基加 N- 前缀并强制省略位次）。"""
    zh_stem = subs[0].get("zh") or ""
    if subs and all((s.get("kind") or "") in N_PREFIX_KINDS for s in subs):  # 整组全为 N-型取代基才走 N-计数前缀；同词干混入 C-型时落入数字通道，N-型成员由 _locant_str 渲染为 N（如 N,N,2-trimethyl）。
        need = any(s.get("paren") for s in subs) or bool(stem and stem[0].isdigit())  # 复合取代基（含 locant 位次/显式 paren）须整体加括号：N-(3-bromophenyl)。
        s_en = _wrap_stem(stem, need)
        s_zh = _wrap_stem(zh_stem, need)
        tokens = _n_prime_tokens(subs, primes)  # 同 N 多取代 → N,N-；跨不同 N → N,N'-（两个甲基挂不同氮时漏撇号会把结构写成另一个分子）
        return _n_prefix_en(len(subs), s_en, tokens, subs), _n_prefix_zh(len(subs), s_zh, tokens, subs)
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit, paren_cf3, stem)


def _n_prime_map(substituents: list) -> dict[int, int]:
    """N-型取代基的 N 原子 → 撇号个数（P-14.5：引用序最前的取代基所在的 N 取未加撇的 N，其余按序加撇）。"""
    groups = _group_by_stem(substituents)
    seen: dict[int, str] = {}
    for stem in _sorted_stems(groups):
        for s in groups[stem]:
            if (s.get("kind") or "") in N_PREFIX_KINDS and s.get("attach_idx") is not None:
                seen.setdefault(s["attach_idx"], stem)
    return {a: i for i, a in enumerate(
        sorted(seen, key=lambda a: (alkyl_alpha_key(seen[a]), a)))}


def _collect_parts(groups: dict[str, list], omit: bool, paren_cf3: bool = False,
                   bracket: bool = False, primes: dict[int, int] | None = None) -> tuple[list[str], list[str]]:
    """汇总所有词干的中英文前缀部件列表；bracket(单碳多不同取代, P-16.5.1.3.1)时第二词干起整体加圆括号(倍增前缀不括入)，词干间无连字符。"""
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for i, stem in enumerate(_sorted_stems(groups)):
        subs = groups[stem]
        if bracket and i >= 1:
            en_parts.append(f"{_mult_en(len(subs))}({stem})")
            zh_parts.append(f"{_mult_zh(len(subs))}({subs[0].get('zh') or ''})")
            continue
        en_p, zh_p = _parts_for_stem(stem, subs, omit, paren_cf3, primes)
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
    en_parts, zh_parts = _collect_parts(groups, omit, paren, bracket, _n_prime_map(substituents))  # P-16.5.1.3.1/.3.2：单碳(meth)母链带 ≥2 个不同简单取代基且位次省略 → 首词干平铺、第二及以后各自括号；单碳链取代基必同处唯一碳，括号式即 locant 省略时的消歧写法。
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
