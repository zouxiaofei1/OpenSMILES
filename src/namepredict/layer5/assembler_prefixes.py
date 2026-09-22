"""L5 取代基前缀分组与双语渲染。"""
from __future__ import annotations

import re

from namepredict.layer1 import fg_registry as _fg_reg
from namepredict.layer4.locant_calc import locant_str_sort
from namepredict.tools.re import alpha_order_key
from namepredict.constants import (
    BIS_EN, BIS_ZH, BRIDGE_DIATOMIC_ZH, BRIDGE_SPLIT_SUFFIX_EN, BRIDGE_SPLIT_SUFFIX_ZH,
    BRIDGE_SUFFIX_EN, BRIDGE_SUFFIX_ZH, CATION_YL_STEMS, DIATOMIC_BRIDGE_YL, MULT_EN, MULT_ZH,
    N_LOCANT_KINDS, N_PREFIX_KINDS, OXO_CENTER_KINDS,
)

def _mult_rows(items: list, key_fn, zh_fn, sort_key=None) -> list[list]:
    """同基分组计数并按 sort_key 排序，返回 [en, zh, 倍数, 成员] 行。"""
    table: dict[str, list] = {}
    for it in items:
        en = key_fn(it)
        row = table.setdefault(en, [en, zh_fn(it), 0, []])
        row[2] += 1
        row[3].append(it)
    rows = list(table.values())
    rows.sort(key=(lambda r: r[0]) if sort_key is None else (lambda r: sort_key(r[0])))
    return rows


def _locant_str(subs: list) -> str:
    """按位次排序拼接成逗号串；N-型取代基渲染为字母位次 N。"""
    tokens = []
    for s in subs:
        if "locant" not in s:
            continue
        kind = s.get("kind") or ""
        tokens.append("N" if kind in N_PREFIX_KINDS else s["locant"])
    return ",".join(str(x) for x in locant_str_sort(tokens))

_DIGIT_RE = re.compile(r"\d")  # 取代基名中的位次数字


def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None,
                      scaffold: str | None = None, has_ene: bool = False) -> bool:
    """判断取代基位次可否省略（环烷烃/苯单取代、酰胺 N- 等情形）。"""
    if kind in N_LOCANT_KINDS:  # 脲/硫脲/胍保留名母体：N 位次须显式写出（1,3-二甲基脲）
        return False
    if kind == "carbamic_acid":  # P-65.2.1.1：N-取代氨基甲酸不带位次（dimethylcarbamic acid）
        return True
    if kind == "radical":  # 自由基母体：连接点隐含 locant 1，单碳链省略位次
        return n_carbons == 1
    if kind == "cation":  # 单核母体阳离子：取代基全挂在同一个原子上，位次恒可省（P-73.1.1）
        return True
    if n_carbons <= 1:  # 单碳母体位次省略；N-/C- 型共存时 C 侧须带位次
        kinds = {(s.get("kind") or "") for s in substituents}
        return not (kinds & N_PREFIX_KINDS and kinds - N_PREFIX_KINDS)
    if (  # 纯烃环单取代位次隐含（环烯须保留位次）
        kind == "alkane" and scaffold in ("carbocycle", "benzene") and not has_ene
    ) and len(substituents) == 1:
        return True
    if kind == "amide":
        return {s.get("kind") for s in substituents} <= N_PREFIX_KINDS
    if kind in ["acyl","ketone","acid"]:
        return False
    if kind in OXO_CENTER_KINDS:  # 含氧酸中心母体：臂挂在中心原子上，母体链位次无意义
        return True
    if any(s.get("paren") or (s.get("en") or "")[:1] == "(" for s in substituents):  # 复合取代基（显式括号）须保留母体 2- 消歧
        return False
    complex_sub = any(  # 自带位次的取代基同须保留，简单 FG 前缀省略
        _DIGIT_RE.search(s.get("en") or "") for s in substituents
    )
    return (  # C2 单取代省略仅对端碳无 H 的母体成立（P-14.3.4.4）
        n_carbons == 2 and len(substituents) == 1
        and kind not in ("alcohol", "amine", "thiol", "sulfonic", "sulfonate",
                         "sulfonamide", "sulfonyl_chloride")  # 端碳已带后缀，C2 取代基位次不可省
        and not complex_sub
    )


_INDICATED_H_LEAD_RE = re.compile(r"^\d+H-")  # 前导指示氢（1H-咪唑-5-基），非取代基位次

def _ring_fused_lead(stem: str) -> bool:
    """词干是否为「环基 + 连接组分」式复合前缀（1H-咪唑-5-基甲基、1,3-苯并二氧杂环戊烯-5-基氧基甲基）。

    此类词干的前导数字属于环系名自身编号（或指示氢），不是取代基位次，
    与母体位次不构成同类位次冲突（P-16.5.1.2）。
    """
    s = stem or ""
    return bool(_INDICATED_H_LEAD_RE.match(s)) or (s[:1].isdigit() and "-yl" in s[:-2])


def _stem_needs_paren(stem: str, subs: list, omit: bool, flat: bool = False) -> bool:
    """加括号：显式标记、前导位次词干、前导立体描述符，或多个 CF3（英文）。"""
    if any(s.get("paren") for s in subs):
        return not _CATION_YL_TAIL_RE.search(stem)  # 单核阳离子基名单词整体（dimethylsulfonio），无需围栏
    if _ring_fused_lead(stem):
        return (not omit) and not flat  # 母体位次省略时无同类位次可混，前缀无需围栏
    if stem and stem[0].isdigit():
        return (not flat) and not (omit and any(s.get("ring_yl") for s in subs))  # 母体位次已省略时环基名无须围栏消歧（P-16.5.1.2）
    if (  # P-16.5.1.4：酰基前缀自带母体氢化物名时须围栏，避免一名两母体；烷氧羰基（…oxycarbonyl）除外
        stem.endswith("carbonyl") and stem != "carbonyl" and "oxy" not in stem
    ):
        return (not omit) and not flat
    return (not omit) and (stem == "trifluoromethyl" or bool(_STEREO_LEAD_ENCLOSE_RE.match(stem)))


def oxo_arm_fence(name: str, sub: dict, mol) -> bool:
    """含氧酸中心母体的臂名围栏：前导立体描述符，或复合基挂在环上。"""
    if (name or "")[:1] == "[":
        return False  # 臂名已自带完整方括号围栏（P-16.5.2），L5 不再二次围栏
    if (name or "")[:1] == "(":
        return True  # 前导圆括号仅为立体描述符（(2S)-…），非完整围栏
    if len(_LOCANT_RUN_RE.findall(_STEREO_LEAD_ENCLOSE_RE.sub("", name or ""))) >= 2:
        return True  # 自带多个位次段的复合臂名须整体围栏（P-16.5.1.3.1）；2,3- 只算一段
    attach = sub.get("attach_idx")
    if not sub.get("paren") or mol is None or attach is None:
        return False
    root = None  # 臂根：claim 原子中与附着原子成键者
    for a in sub.get("atoms") or ():
        try:
            if mol.GetBondBetweenAtoms(int(a), int(attach)) is not None:
                root = int(a)
                break
        except (ValueError, TypeError):
            return False
    if root is None:
        return False
    return mol.GetAtomWithIdx(root).IsInRing()


def _sub_root_in_ring(sub: dict, mol) -> bool:
    """取代基根原子（与附着原子成键者）是否在环上（P-16.5.1.2）。"""
    attach = sub.get("attach_idx")
    if attach is None:
        return False
    for a in sub.get("atoms") or ():
        try:
            if mol.GetBondBetweenAtoms(int(a), int(attach)) is not None:
                return mol.GetAtomWithIdx(int(a)).IsInRing()
        except (ValueError, TypeError):
            return False
    return False


def _enclose(s: str) -> str:
    """给名称加围栏：已含圆括号时改用方括号（嵌套规则），否则加圆括号。"""
    return f"[{s}]" if "(" in s else f"({s})"


def _wrap_stem(stem: str, need: bool) -> str:
    """按 need 给词干加括号（词干含括号时改用方括号）。"""
    return _enclose(stem) if need else stem


def _place(mult: str, s: str, subs: list, omit: bool) -> str:
    """拼数量前缀与词干体：位次省略时直接相接，否则位次串以连字符前置。"""
    if omit:
        if mult and s[:1].isdigit():  # 倍数前缀不得与位次数字直连（P-16.3.2：bis(3-…)，非 bis3-…）
            return f"{mult}({s})"
        return f"{mult}{s}"
    return f"{_locant_str(subs)}-{mult}{s}"

# 取代基名以带位次的立体描述符开头：(1Z)-、(2R,4R)-。与 tools.re._STEREO_LEAD_STRIP_RE 不同——
# 那份供 P-14.5 排序剥除、允许无位次，本份供围栏判定、要求每位次带 token，勿互换。
_STEREO_LEAD_ENCLOSE_RE = re.compile(r"\(\d+[RSEZ](?:,\d+[RSEZ])*\)-")

_BRIDGE_SELF_FENCE = ("sulfanyl", "sulfinyl", "硫基", "亚磺酰基")  # 前端自带方括号时并入同一围栏的桥后缀（P-16.5.1.3）
_FRONT_TAILS = ("yl", "ylidene", "sulfanyl", "amino")  # 可作桥前端的词尾：-yl/-ylidene 自由价基，或本身即复合桥前端（…amino）
_MULT_WRAP_RE = re.compile(r"-\d+-yl$|oyloxy$")  # 须整体加括号的复合词干：位次链基与酰氧基


_CHAIN_STEM = (r"(?:meth|eth|prop|but|pent|hex|hept|oct|non|dec|undec|dodec|tridec|tetradec|"
               r"pentadec|hexadec|heptadec|octadec|nonadec|eicos)")
_SIMPLE_CHAIN_YL_RE = re.compile(r"^(?:\d+-)?" + _CHAIN_STEM + r"a?n-\d+-yl$")  # 无取代直链 -yl 与 O/S 桥融合平铺，不拆围栏。
_SUBST_CHAIN_YL_RE = re.compile(_CHAIN_STEM + r"a?n-\d+-yl$")  # 带取代基的直链 -yl 仍与桥融合平铺。
_TERMINAL_CHAIN_YL_RE = re.compile(_CHAIN_STEM + r"yl$")  # 自由价在端碳的直链基与桥融合平铺。
_BENZYL_TAIL_RE = re.compile(r"[)\]]methyl$")  # 苄基型前端：桥后缀直接缀在甲基上，不拆。
_LOCANT_RE = re.compile(r"(?:^|[-,\[])\d")  # 位次数字：行首或 -,\[ 之后（立体描述符内的数字不算）
_ACYL_FRONT_RE = re.compile(r"(?:oyl|carbonyl)$")  # 酰基前端词尾（乙酰氧/酰胺融合用）
_LOCANT_SUBST_TAIL_RE = re.compile(r"[\d\]]-[^()]*yl$")  # 括号外仍带位次取代基的端基（1-(…)-4-methylsulfanylbutyl）
_CATION_YL_TAIL_RE = re.compile(r"(onio|onium|inium)$")  # 单核阳离子去氢基名词尾（dimethylsulfonio）
_MERGE_TAIL_RE = re.compile(r"\][a-z]*yl$")  # 端部为「方括号组 + 直链 -yl」：桥后缀可并入同一围栏


_RING_STEM_RE = re.compile(  # 常见环系基名词干（判定复合前缀用）
    r"cyclo|benzen|naphthalen|anthracen|phenanthren|pyridin|pyrimidin|pyridazin|pyrazin|triazin|"
    r"pyrrol|imidazol|pyrazol|triazol|tetrazol|oxazol|thiazol|oxan|thian|oxolan|thiolan|furan|thiophen|"
    r"indol|indazol|indolin|quinolin|isoquinolin|purin|pteridin|chromen|chromanon|piperidin|piperazin|"
    r"morpholin|azepan|azetidin|aziridin|oxepan|thiepan|picen|gonan|androstan|estran|pregnan|cholestan"
)


_FRONT_LOCANT_NUM_RE = re.compile(r"(?:^|[-,\[])(\d+)")  # 前端名内位次数字（不含立体描述符内的数字）


def _locant_conflict(base: str) -> bool:
    """同一位次数字在前端名内出现二次：指向不同结构要素，须围栏（P-16.5.1.2）。"""
    nums = _FRONT_LOCANT_NUM_RE.findall(base)
    return len(nums) != len(set(nums))


def _composite_front(base: str) -> bool:
    """桥前端是否为复合前缀（P-16.5.1.1）：环系基上另挂烃基取代基者须围栏。"""
    hit = _RING_STEM_RE.search(base)
    return bool(hit) and hit.start() > 0 and bool(_LOCANT_RE.search(base))


def _front_needs_enclosure(base: str, suf: str) -> bool:
    """O/S/N 桥前端是否为自带围栏的复合取代基。"""
    if base.endswith(("sulfanyl", "amino")):  # 前端自身即复合桥名（…aminooxy）：后端另起一重前缀（P-63.2.2.1）
        return True
    if base.endswith(("sulfonyl", "sulfinyl")):  # 磺酰基前端围栏由 L3 定形，此处不拆。
        return False
    if _BENZYL_TAIL_RE.search(base):  # 苄基型前端（…yl]methylsulfanyl）：桥后缀直接缀在甲基上，不拆。
        return False
    if base.endswith("ylidene") and "[" not in base:  # 简单 (X)亚基前端与桥融合（…methylideneamino）
        return False
    if suf == "oxy" and base.endswith("carbonyl"):  # P-63.2.2.1：酰基前端与氧桥融合（…carbonyloxy）
        return False
    if suf == "amino" and _ACYL_FRONT_RE.search(base):  # P-63.2.2.1.2：amino 桥酰基前端按取代式融合（…oylamino）
        return False
    if  re.match(r"^\(\d+[RrSs]", base) and not base.endswith("oyl"):  # 手性自由价碳前端须括起，酰基前端按 …oyloxy 融合
        return True
    if "[" in base:  # 方括号前端：amino 桥按位次前缀细分，其余不拆
        return True if suf != "amino" else bool(
            re.search(r"\]-?\d", base)                            # 括号后接数字位次前缀
            or re.search(r"-\d+-\[", base[: base.find("[") + 1]))  # 括号前已有数字位次前缀：3-oxo-3-[X]propyl
    if "(" in base:  # 前端自带括号：端碳自由价链基平铺，同一位次指向不同要素者须围栏（P-16.5.1.2）。
        flat = _TERMINAL_CHAIN_YL_RE.search(base) and not _locant_conflict(base)
        return (not flat) or bool(_LOCANT_SUBST_TAIL_RE.search(base))
    if _composite_front(base):  # 复合前缀（自带位次的环基/链基）作桥前端须围栏（P-16.5.1.1）
        return True
    if suf in ("oxy", "sulfanyl", *DIATOMIC_BRIDGE_YL) and re.search(r"\d", base) and base.endswith("phenyl"):  #
        return True
    if suf == "amino":  # P-63.2.2.1.2：amino 桥按取代式融合，不拆。
        return False
    return bool(re.search(r"-\d+-yl$", base)) and not _SUBST_CHAIN_YL_RE.search(base)


def _split_bridge_suffix(stem: str) -> tuple[str, str] | None:
    """拆 O/S/N 桥平铺式为 (前端, 桥后缀)（P-63.2.2.1.1）。"""
    for suf in BRIDGE_SPLIT_SUFFIX_EN:
        if not stem.endswith(suf):
            continue
        base = stem[: -len(suf)]
        if not base.endswith(_FRONT_TAILS) or not _front_needs_enclosure(base, suf):
            continue
        return base, suf
    return None


def _sbridge_flat_stem(stem: str) -> bool:
    """磺酰基/亚磺酰基桥 + 直链 -yl 前端英文平铺（P-63.2.1）。"""
    return any(stem.endswith(suf) and _SIMPLE_CHAIN_YL_RE.match(stem[: -len(suf)])
               for suf in ("sulfonyl", "sulfinyl"))


def _bridge_body(base: str, suf: str, merge: bool = False) -> str:
    """O/S/N 桥平铺主体：前端围栏 + 桥后缀留外（P-63.2.2.1.2）；merge 时整段自闭合围栏。"""
    body = f"{_enclose(base)}{suf}"
    if merge:  # 该前缀后仍接其他前缀：括界须自行闭合（P-16.5.1.3.1）
        return _enclose(f"{base}{suf}") if _MERGE_TAIL_RE.search(base) else f"[{body}]"
    if suf in DIATOMIC_BRIDGE_YL + BRIDGE_DIATOMIC_ZH:  # 双原子桥：前端已围栏，整段再括一层
        return f"[{body}]"
    if suf in ("amino", "氨基"):  # 前端围栏 + 氨基：整段再括一层（P-16.5.4.1.5 括号种类升级）
        return f"[{body}]"
    return body


def _prefix_one_en(stem: str, subs: list, omit: bool, tail_sep: bool = False,
                   flat: bool = False) -> str:
    """拼单个英文前缀：数量 + 词干（可省略位次时省略 locant）。"""
    mult = _mult_of("en", stem, subs, len(subs))
    need = (_stem_needs_paren(stem, subs, omit, flat)
            and not _is_bare(subs) and not _sbridge_flat_stem(stem))
    if mult and mult == MULT_EN.get(len(subs), "") and _MULT_WRAP_RE.search(stem):  # P-16.3.2：复合取代基的倍数前缀须加括号
        need = True
        mult = "bis" if stem[:1] in "aeiou" else "di"
    if len(subs) > 1 and stem in CATION_YL_STEMS:  # P-16.3.2：阳离子去氢前缀属复合前缀，倍数用 bis(...)：bis(azaniumyl)
        need, mult = True, BIS_EN.get(len(subs), mult)
    if need:
        sp = _split_bridge_suffix(stem)
        if sp is not None:  # O/S/N 桥平铺式：桥后缀留括号外（P-63.2.2.1.1）。
            merge = tail_sep and sp[1] in _BRIDGE_SELF_FENCE and "[" in sp[0]
            body = _bridge_body(*sp, merge=merge)
            if len(subs) > 1 and mult in (MULT_EN.get(len(subs)), BIS_EN.get(len(subs))):
                return _place(BIS_EN.get(len(subs), mult), _enclose(body), subs, omit)  # 桥后缀留在围栏外：倍数组须整体围栏（P-16.3.2）
            return _place(mult, body, subs, omit)
    return _place(mult, _wrap_stem(stem, need), subs, omit)

_COMPLEX_MULT_LANG = {"en": ("carboxy", BIS_EN, MULT_EN), "zh": ("羧", BIS_ZH, MULT_ZH)}


def _mult_of(lang: str, stem: str, subs: list, n: int) -> str:
    """数量前缀（P-16.3.2）：复合组分用 bis/tris。"""
    sentinel, bis, mult = _COMPLEX_MULT_LANG[lang]
    if sentinel in stem or any(s.get("paren") for s in subs):
        return bis.get(n, "")
    return mult.get(n, "")


def _split_bridge_suffix_zh(zh_stem: str, en_stem: str) -> tuple[str, str] | None:
    """中文侧 O/S/N 桥平铺式拆分：判据与英文侧同步。"""
    sp = _split_bridge_suffix(en_stem)
    if sp is None:
        return None
    for suf in BRIDGE_SPLIT_SUFFIX_ZH:
        if zh_stem.endswith(suf):
            base = zh_stem[: -len(suf)]
            if base.endswith("基"):
                return base, suf
            if base.endswith("-") or base.endswith("氨"):  # 桥融合时「基」被氧基顶掉（…氨氧基），拆时补回。
                return f"{base}基", suf
    return _zh_front_ji(zh_stem, sp)


def _zh_front_ji(zh_stem: str, sp: tuple[str, str]) -> tuple[str, str] | None:
    """复合前端的「基」被上游切掉：前端自带围栏时在拆分点补回（酰基前端走融合式除外）。"""
    if _ACYL_FRONT_RE.search(sp[0]) or not _front_needs_enclosure(*sp):
        return None
    for suf in BRIDGE_SPLIT_SUFFIX_ZH:
        base = zh_stem[: -len(suf)]
        if zh_stem.endswith(suf) and base and not base.endswith("基"):
            return f"{base}基", suf
    return None


def _prefix_one_zh(zh_stem: str, subs: list, omit: bool,
                   en_stem: str = "", tail_sep: bool = False, flat: bool = False) -> str:
    """拼单个中文前缀：数量 + 词干（CF3 特例：简单氟代甲基不加括号）。"""
    mult = _mult_of("zh", zh_stem, subs, len(subs))
    en = subs[0].get("en") or ""
    omit_ring_yl = omit and any(s.get("ring_yl") for s in subs)  # 母体位次已省时环基词干免围栏，与英文侧同判（P-16.5.1.2）
    need = (any(s.get("paren") for s in subs)
            or bool(en[:1].isdigit() and not flat and not omit_ring_yl)) and not _is_bare(subs)
    if not omit and _STEREO_LEAD_ENCLOSE_RE.match(en):  # 前导立体描述符 + 位次须整体围栏
        need = True
    if mult and mult == MULT_ZH.get(len(subs), "") and en_stem and _MULT_WRAP_RE.search(en_stem):  # 与英文侧同判：复合取代基用 双(...)
        need = True
        mult = BIS_ZH.get(len(subs), mult)
    if len(subs) > 1 and en_stem in CATION_YL_STEMS:  # 与英文侧同判：阳离子去氢前缀用 双(铵基)
        need, mult = True, BIS_ZH.get(len(subs), mult)
    if need:
        sp = _split_bridge_suffix_zh(zh_stem, en_stem)
        if sp is not None:  # 与英文侧同形：桥后缀留括号外（P-63.2.2.1.1）
            en_sp = _split_bridge_suffix(en_stem)
            merge = tail_sep and bool(en_sp) and en_sp[1] in _BRIDGE_SELF_FENCE and "[" in en_sp[0]
            return _place(mult, _bridge_body(*sp, merge=merge), subs, omit)
    return _place(mult, _wrap_stem(zh_stem, need), subs, omit)


def _n_prefix(lang: str, n: int, stem: str, tokens: list[str], subs: list) -> str:
    """N- 前缀：N-甲基 / N,N-二甲基（N,N'-双）。"""
    if n == 1:
        return f"{tokens[0]}-{stem}"
    return f"{','.join(tokens)}-{_mult_of(lang, stem, subs, n)}{stem}"


def _is_bare(subs: list) -> bool:
    """该组前缀是否免去二次围栏（单核母体阳离子的臂名前置于母体名，无位次可混，P-73.6）。"""
    return any(s.get("bare") for s in subs)


def _strip_nested_fence(name: str) -> str:
    """由内向外剥掉嵌套围栏，只留顶层正文（用于数顶层位次段）。"""
    while True:
        stripped = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", "", name)
        if stripped == name:
            return name
        name = stripped


def cation_arm_bare(name: str) -> bool:
    """阳离子母体臂名是否免围栏：前导立体描述符或自带多个位次段者须围栏。"""
    if _STEREO_LEAD_ENCLOSE_RE.match(name or ""):
        return False
    top = _strip_nested_fence(name or "")
    return len(_LOCANT_RUN_RE.findall(top)) < 2  # 2,3- 只算一段


def _parts_for_stem(stem: str, subs: list, omit: bool, primes: dict[int, int] | None = None,
                    tail_sep: bool = False, flat: bool = False) -> tuple[str, str]:
    """按词干生成中英文前缀（N- 类取代基加 N- 前缀并强制省略位次）。"""
    zh_stem = subs[0].get("zh") or ""
    if subs and all((s.get("kind") or "") in N_PREFIX_KINDS for s in subs):  # 整组全为 N-型才走 N-计数前缀，混入 C-型时走数字通道
        need = any(s.get("paren") for s in subs) or bool(stem and stem[0].isdigit())  # 复合取代基须整体加括号
        s_en = _wrap_stem(stem, need)
        s_zh = _wrap_stem(zh_stem, need)
        tokens = sorted(("N" + "'" * ((primes or {}).get(s.get("attach_idx"), 0))) for s in subs)  # 同 N 用 N,N-；跨不同 N 用 N,N'-
        return (_n_prefix("en", len(subs), s_en, tokens, subs),
                _n_prefix("zh", len(subs), s_zh, tokens, subs))
    return (_prefix_one_en(stem, subs, omit, tail_sep, flat),
            _prefix_one_zh(zh_stem, subs, omit, stem, tail_sep, flat))


def _n_prime_map(groups: dict[str, list], stems: list[str]) -> dict[int, int]:
    """N-型取代基的 N 原子 → 撇号个数（P-14.5 引用序）。"""
    seen: dict[int, str] = {}
    for stem in stems:
        for s in groups[stem]:
            if (s.get("kind") or "") in N_PREFIX_KINDS and s.get("attach_idx") is not None:
                seen.setdefault(s["attach_idx"], stem)
    return {a: i for i, a in enumerate(
        sorted(seen, key=lambda a: (alpha_order_key(seen[a]), a)))}


def _arm_tail(name: str, n: int, bare: bool, mult: dict, bis: dict) -> str:
    """括号式后继臂（P-16.5.1.3.1）：原式「倍数(名)」；阳离子臂按 P-73.1.1 改直连式。"""
    if not bare:
        return f"{mult.get(n, '')}({name})"
    if name[:1].isdigit() or "(" in name or "[" in name:  # 带位次/围栏的复合臂：倍数前缀须整体成段
        return f"{bis.get(n, '')}({name})" if n > 1 else _enclose(name)
    return f"{mult.get(n, '')}{name}" if n > 1 else name


def _collect_parts(groups: dict[str, list], stems: list[str], omit: bool,
                   bracket: bool = False, primes: dict[int, int] | None = None,
                   sep: str = "-", flat: bool = False,
                   bare: bool = False) -> tuple[list[str], list[str]]:
    """汇总所有词干的中英文前缀部件列表（P-16.5.1.3.1 括号式）。"""
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for i, stem in enumerate(stems):
        subs = groups[stem]
        tail_sep = sep == "-" and i < len(stems) - 1  # 该前缀后仍接别的前缀（括界须自行闭合）
        if bracket and i >= 1:
            en_parts.append(_arm_tail(stem, len(subs), bare, MULT_EN, BIS_EN))
            zh_parts.append(_arm_tail(subs[0].get("zh") or "", len(subs), bare, MULT_ZH, BIS_ZH))
            continue
        if bracket and not _LOCANT_RE.search(stem):  # P-16.5.1.3.1：首个引用的取代基从不加围栏（自带位次者除外）
            subs = [{**s, "paren": False} for s in subs]
        en_p, zh_p = _parts_for_stem(stem, subs, omit, primes, tail_sep, flat)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts


_LOCANT_RUN_RE = re.compile(r"\d+(?:,\d+)*[a-z]*-")  # 名内位次段（2,3- 算一段），与 tools.re.SUB_LOCANT_RE 不同：无段首锚点


def _cation_arm_hyphen(groups: dict[str, list]) -> bool:
    """阳离子臂名间是否以连字符分段（P-73.1.1）：三臂以上，或存在复合臂时须显式分隔。"""
    if len(groups) >= 3:
        return True
    for name, subs in groups.items():
        if "(" in name or "[" in name:  # 自带围栏的复合臂
            return True
        if name[:1].isdigit() and (len(_LOCANT_RUN_RE.findall(_strip_nested_fence(name))) >= 2
                                   or len(subs) > 1):  # 多位次臂名，或带位次的倍增臂
            return True
    return False


def _groups_simple(groups: dict[str, list]) -> bool:
    """全部词干无前导位次且非 N- 类才适用括号式（前导位次词干已由位次连字符式消歧）。"""
    for subs in groups.values():
        for s in subs:
            if (s.get("en") or "")[:1].isdigit() or (s.get("kind") or "") in N_PREFIX_KINDS:
                return False
    return True


_O_SIDE_ARM_FENCE_KINDS = frozenset({"sulfonate", "phosphate"})  # 组装侧不判 O-侧臂围栏的两类中心母体


def _o_side_arm_fence(name: str, sub: dict) -> bool:
    """O-侧臂围栏（P-16.5.1.3.1）：自带多位次的复合臂名须整体括起。"""
    if not sub.get("paren") or (name or "")[:1] in "[(" or re.search(r"[()\[\]]", name or ""):
        return False  # 简单臂名与已自带围栏/含括号的臂名（立体描述符、复合前缀）不加
    return len(_LOCANT_RUN_RE.findall(name)) >= 2


def _fence_o_side_arms(subs: list, kind: str | None, mol) -> None:
    """酯路径 O-侧臂名就地围栏：assembler 侧对这两类 kind 不判臂围栏。"""
    if kind not in _O_SIDE_ARM_FENCE_KINDS or mol is None:
        return
    for s in subs:
        if not s.get("o_side") or s.get("arm_fenced"):
            continue
        s["arm_fenced"] = True  # 同一 numbered 会被多次组装，标记防二次围栏
        name = s.get("en") or ""
        attach = s.get("attach_idx")
        if not name or attach is None:
            continue
        if mol.GetAtomWithIdx(int(attach)).GetAtomicNum() == 16:  # 硫代酯 S-侧臂：assembler 侧已围栏
            continue
        if _o_side_arm_fence(name, s):
            s["en"] = _enclose(name)


def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None,
                  scaffold: str | None = None, has_ene: bool = False,
                  mol=None, locant_kind: str | None = None) -> tuple[str, str]:
    """组合完整取代基前缀：滤 O 侧、判 omit、按词干分组拼接。"""
    _fence_o_side_arms(substituents, kind, mol)
    substituents = [s for s in substituents if not s.get("o_side")]  # ester 的 O 侧臂由 join_kind_name 消费
    if not substituents:
        return "", ""
    if mol is not None:  # 环基取代基名自带位次体系，标记后供 _stem_needs_paren 免围栏
        substituents = [{**s, "ring_yl": _sub_root_in_ring(s, mol)} for s in substituents]
    omit = _omit_sub_locants(n_carbons, substituents, locant_kind or kind, scaffold, has_ene)
    bare = kind == "cation"  # 单核母体阳离子：臂名直接前置于阳离子名，无母体位次可混，不再二次围栏
    if bare:
        substituents = [{**s, "bare": cation_arm_bare(s.get("en") or "")} for s in substituents]
    flat = kind in OXO_CENTER_KINDS  # 中心母体：臂名围栏改由环基判据决定，去多余位次围栏
    if flat:
        substituents = [{**s, "paren": oxo_arm_fence(s.get("en") or "", s, mol)} for s in substituents]
    rows = _mult_rows(substituents, lambda s: s.get("en") or "", lambda s: s.get("zh") or "",
                      alpha_order_key)  # P-14.5：全部前缀按字母数字序引用，非仅烷基
    groups = {en: members for en, _, _, members in rows}
    stems = [en for en, _, _, _ in rows if en]
    bracket = bool(omit) and len(groups) >= 2 and (bare or _groups_simple(groups))  # P-16.5.1.3.1：位次省略的单核母体，首基平铺、余基括起
    arm_hyphen = bool(bare and bracket and _cation_arm_hyphen(groups))  # 阳离子臂名整体前置于阳离子名，复合臂间以连字符分段（P-73.1.1）
    sep = "-" if (arm_hyphen or not bracket) else ""
    en_parts, zh_parts = _collect_parts(groups, stems, omit, bracket, _n_prime_map(groups, stems), sep, flat, arm_hyphen)  # P-16.5.1.3.1/.3.2：单碳链多不同取代基 → 首平铺，余加括号
    return sep.join(en_parts), sep.join(zh_parts)

def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    """从 numbered 提取上下文并委托 _build_prefix。"""
    parent = numbered.get("parent") or {}
    has_ene = bool(parent.get("double_bond") or parent.get("double_bonds"))
    if parent.get("subs_consumed"):  # 取代基已并入 C1 保留名（carbamoyl 等），不再另加前缀
        return "", ""
    if kind == "radical" and parent.get("radical_anchor_element"):  # 杂原子锚点自由基：烷基取代基已并入组装名（ethyloxy），不再加前缀。
        return "", ""
    return _build_prefix(numbered.get("substituents") or [], n, kind,
                         parent.get("scaffold_id"), has_ene, parent.get("mol"),
                         locant_kind=parent.get("locant_kind"))
