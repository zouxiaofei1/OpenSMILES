"""L5 名称组装主入口：由链引擎取名后按 kind 拼接前缀、立体（E/Z、R/S）与盐类后缀。"""
from __future__ import annotations
from dataclasses import replace
import re
from namepredict.constants import (
    ALKOXY_YLOXY_EN, ALKOXY_YLOXY_ZH, AMIDO_RETAINED, AZANE_PAREN_SUF, BIS_EN, BRIDGE_YL_SUFFIX,
    BRIDGE_ZH_YL_SUFFIX,
    EXO_RING_SUF, MONONUCLEAR_BRIDGE, MONONUCLEAR_YL, MONONUCLEAR_ZERO_YL, MULT_EN,
    MULT_ZH, PHOSPHORYL_STEMS, zh_bridge_root,
)
from namepredict.layer5.chain_engine import _ACYL_HALIDE_BY_HAL, _KIND_TABLE, _alkane_names, _chain_names
from namepredict.layer5.stems import join_anion_names, join_metal_salt_names
from namepredict.layer5.assembler_prefixes import _SIMPLE_CHAIN_YL_RE, _prefix_for
from namepredict.layer5.stereo import _split_stereo_lead as _stereo_lead
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult:
    """构造失败 NameResult（success=False，meta 供诊断）。"""
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    """构造成功 NameResult（success=True，记录耗时与来源）。"""
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

# 稠环/杂环词干（en_stem, zh_stem）：由 L2 注入 parent 的 stem_en/stem_zh 派生。IUPAC 词干 = 母体名去尾部 e（benzene 例外，保留完整名）；aromatic 恒 True（保留母体均芳香）。
def _ring_stem(numbered: dict) -> tuple[str, str] | None:
    """保留 scaffold 的完整 IUPAC 词干（用于 -ol/-diol/-amine 等 FG 后缀拼接）；结尾 'e' 的省略交给 chain_engine._elide_parent_e 按后缀首字母判断（P-60.2(a)），避免 oxolane-3,4-diol 被错拼成 oxolan-3,4-diol。"""
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None


# 稠环/杂环完整 base 名（-carboxylic acid 用完整词干，如 naphthalene-1-carboxylic acid）。
def _ring_base(numbered: dict) -> tuple[str, str] | None:
    """保留 scaffold 的完整母体名。"""
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None


def _scaffold_id(numbered: dict) -> str | None:
    """取母体的 scaffold_id（carbocycle/benzene/稠环…）。"""
    return (numbered.get("parent") or {}).get("scaffold_id")


def _ring_carbocycle_stem(n: int, numbered: dict) -> tuple[str | None, str | None, bool]:
    """单环 carbocycle 母体词干（含环内烯不饱和度）：返回 (en_stem, zh_stem, has_unsat)；en 单烯 1 位省略（P-31.1.2）、多烯/非 1 位显式，中文位次恒显式，has_unsat=True 时调用方须显式给出 FG 位次（P-65/P-66 中烯使环编号不再唯一），en/zh 为 None 表示不支持。"""
    base = _alkane_names(n)
    if not base:
        return None, None, False
    enes = numbered.get("ene_locants")
    if not enes:
        return f"cyclo{base[0]}", f"环{base[1]}", False
    en_core = base[0][:-3] if base[0].endswith("ane") else base[0]   # hexane → hex
    zh_core = base[1][:-1] if base[1].endswith("烷") else base[1]    # 己烷 → 己
    if len(enes) >= 2:
        loc = ",".join(str(x) for x in enes)
        m_en, m_zh = MULT_EN.get(len(enes)), MULT_ZH.get(len(enes))
        if not m_en or not m_zh:
            return None, None, True
        return (f"cyclo{en_core}a-{loc}-{m_en}ene", f"环{zh_core}-{loc}-{m_zh}烯", True)
    if enes[0] == 1:
        return f"cyclo{en_core}ene", f"环{zh_core}-1-烯", True
    return f"cyclo{en_core}-{enes[0]}-ene", f"环{zh_core}-{enes[0]}-烯", True


def _ring_extra_prefix_located(numbered: dict) -> bool:
    """单环环烷上是否另带被编号前缀（oxo/烷基/卤素/羟基…）；exocyclic 主基后缀锚定 locant 1 后，任何其它环位取代（O 侧酯烷基、N 端胺取代除外）必带前缀位次，此时后缀 locant 不可省略（P-66.6.1：4-formylcyclohexane-1-carboxylic acid vs 无取代的 cyclohexanecarbaldehyde）。"""
    chain = set((numbered.get("parent") or {}).get("chain") or [])
    for s in numbered.get("substituents") or []:
        if s.get("o_side"):
            continue
        a = s.get("attach_idx")
        if a is not None and a in chain:
            return True
    return False
def _exocyclic_ring_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """环外主基（羧基/醛/酯/酰胺/腈/酰基头）系统名：环母体词干 + `EXO_RING_SUF` 后缀组装（P-65.1.7.2 酰基头 / P-65.2.2 多羧酸 / P-66.6.1.1.3 环醛）；苯单取代走 chain_engine 保留名（benzoic acid/benzaldehyde/benzonitrile…）时返回 None 回落。"""
    parent = numbered.get("parent") or {}
    facts = parent.get("principal_expression_facts")
    sid = parent.get("scaffold_id")
    spec = EXO_RING_SUF.get(facts.group_class.value) if facts else None
    if spec is None or facts.relation.value != "exocyclic":
        return None
    singular, plural, plural_needs_loc = spec
    mult = facts.multiplicity
    if mult == 1:
        if sid == "benzene":  # 苯单取代保留名由 chain_engine variant 承担，此处不拼 base- 通用名。
            return None
        suf_en, suf_zh = singular
    else:
        if plural is None:  # 酯/酰胺/腈/酰基头无多取代系统名。
            return None
        m_en, m_zh = MULT_EN.get(mult), MULT_ZH.get(mult)
        if not m_en or not m_zh:
            return None
        suf_en, suf_zh = f"{m_en}{plural[0]}", f"{m_zh}{plural[1]}"
    rec = next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == facts.group_class.value), None)
    locs = rec.get("locants") if rec else None
    loc = ",".join(str(x) for x in locs) if locs else None
    if mult > 1 and plural_needs_loc and not loc:  # 多取代位次必带（P-65.2.2），缺失即放弃。
        return None
    if sid == "carbocycle" and not parent.get("fused_tree"):  # 未注册全碳稠环(carbocycle 兜底 + fused_tree)走下方 base 分支, 不作单环环烷烃命名。
        en_ring, zh_ring, has_unsat = _ring_carbocycle_stem(n, numbered)
        if en_ring is None:
            return None
        use_loc = bool(loc) and (mult > 1 or has_unsat or _ring_extra_prefix_located(numbered))  # 单取代：环烯使编号不再唯一、或环上另带前缀取代（必带位次）时主基 locant 1 不可省略（P-65.2.2.1/P-66.6.1）；多取代恒带位次。
        if use_loc:
            return (f"{en_ring}-{loc}-{suf_en}", f"{zh_ring}-{loc}-{suf_zh}")
        return (f"{en_ring}{suf_en}", f"{zh_ring}{suf_zh}")
    base = _ring_base(numbered)
    if base:  # 主基位次取 L4 已算出的 FG locant；无则省略（1 位隐含）。
        if loc:
            return (f"{base[0]}-{loc}-{suf_en}", f"{base[1]}-{loc}-{suf_zh}")
        return (f"{base[0]}{suf_en}", f"{base[1]}{suf_zh}")
    return None


def _bracket_bridge_suffix(en: str, zh: str) -> tuple[str, str]:
    """复合组分加方括号；-yl 型烷氧/硫基把桥后缀挪到括号外（[(2R,…)环己基]氧基，P-16.5.2 嵌套标记）。"""
    for yl_suf, bridge in BRIDGE_YL_SUFFIX:
        if en.endswith(yl_suf):
            en = f"[{en[: -len(bridge)]}]{bridge}"
            break
    else:
        en = f"[{en}]"
    for yl_suf, bridge in BRIDGE_ZH_YL_SUFFIX:
        if zh.endswith(yl_suf):
            zh = f"[{zh[: -len(bridge)]}]{bridge}"
            break
    else:
        zh = f"[{zh}]"
    return en, zh


def _oxido_arm(s: dict, mol) -> tuple[str, str] | None:
    """P 上氧负离子臂（–O⁻）的取代基名 oxido/氧化（P-72.6.2；表 4-2「羟基(氧负离子基)膦酰基 hydroxyoxidophosphoryl」），非阴离子单氧臂返回 None。"""
    atoms = s.get("atoms") or []
    if mol is None or len(atoms) != 1:
        return None
    a = mol.GetAtomWithIdx(int(atoms[0]))
    if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
        return ("oxido", "氧化")
    return None


def _phosphoryl_sub_names(subs: list[dict], stem_en: str, stem_zh: str, mol=None) -> tuple[str, str] | None:
    """P 酰基前缀的取代基拼接（P-67.1.4.1.1.5）：取代基按字母序接到 phosphoryl；全为简单基时首基平铺、其余括起（hydroxy(methyl)phosphoryl），同基倍增用 di-/tri-（dimethoxyphosphoryl）；含复合组分时逐组分以连字符分隔、需围栏者加方括号（P-16.5.2 嵌套标记）。"""
    from namepredict.tools.re import alkyl_alpha_key

    groups: dict[str, list] = {}
    for s in subs:
        en, zh = (s.get("en") or "").strip(), (s.get("zh") or "").strip()
        if not en or not zh:
            return None
        oxido = _oxido_arm(s, mol)
        if oxido is not None:  # 酸式 H 已被夺去的 O⁻ 臂：hydroxy → oxido
            en, zh = oxido
        row = groups.setdefault(en, [en, zh, 0])
        row[2] += 1
    rows = sorted(groups.values(), key=lambda t: alkyl_alpha_key(t[0]))
    if len(rows) == 1 and rows[0][2] > 1:  # 同基倍增：dimethoxyphosphoryl（P-16.3.2 简单基用 di-，不逐基加括号）
        en, zh, m = rows[0]
        m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
        if not m_en or not m_zh:
            return None
        return f"{m_en}{en}{stem_en}", f"{m_zh}{zh}{stem_zh}基"
    compound = any("(" in en or "[" in en for en, _, _ in rows)  # 存在自身带括号/方括号的复合组分：改逐组分连字符分隔 + 方括号围栏（P-16.5.2 嵌套），否则按简单组分平铺/括号
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for i, (en, zh, m) in enumerate(rows):
        if m > 1:
            m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
            if not m_en or not m_zh:
                return None
            en, zh = f"{m_en}{en}", f"{m_zh}{zh}"
        if not compound:  # 简单组分：首基平铺、其余括起
            en_parts.append(en if i == 0 else f"({en})")
            zh_parts.append(zh if i == 0 else f"({zh})")
            continue
        en_need = "(" in en and (i > 0 or en.startswith("("))  # 首组分仅在前导括号（立体描述符）时才须围栏；后续组分带括号即围栏
        zh_need = "(" in zh and (i > 0 or zh.startswith("("))
        if en_need or zh_need:
            br_en, br_zh = _bracket_bridge_suffix(en, zh)
            en, zh = (br_en if en_need else en), (br_zh if zh_need else zh)
        en_parts.append(en)
        zh_parts.append(zh if i == 0 else ("-" if zh.startswith("[") else "") + zh)
    joiner = "-" if compound else ""
    return joiner.join(en_parts) + stem_en, "".join(zh_parts) + stem_zh + "基"


def _azane_acyl_stereo_lead(en: str) -> bool:
    """单 N-酰基残基名是否带前导立体描述符且为酰基词干（azane 方法 2 需括起 acyl 再缀 amino）。"""
    if not en.startswith("("):
        return False
    tag, stem = _stereo_lead(en)
    return bool(tag) and (stem.endswith("oyl") or "carbonyl" in stem)

def _azane_sub_needs_paren(a: dict) -> bool:
    """azane 单取代基是否需整体括起再缀 amino（P-16.5.1.1）：仅限复合标记 + 芳香酰基/acetyl 系内层组。"""
    if not a.get("paren"):  # 简单取代基（methyl/chloro…）无歧义，平铺。
        return False
    return (a.get("en") or "").endswith(AZANE_PAREN_SUF)


def _retained_alkoxy(en: str, zh: str) -> tuple[str, str]:
    """O 锚点 -yloxy 尾部收拢为 IUPAC 保留烷氧基（未命中原样返回）。"""
    for suf_en, kept_en in ALKOXY_YLOXY_EN:
        if en.endswith(suf_en):
            en = en[: -len(suf_en)] + kept_en
            break
    for suf_zh, kept_zh in ALKOXY_YLOXY_ZH:
        if zh.endswith(suf_zh):
            zh = zh[: -len(suf_zh)] + kept_zh
            break
    return (en, zh)


def _bridge_enclosed_names(a: dict, stem_en: str, stem_zh: str) -> tuple[str, str] | None:
    """S 桥前端为复合取代基时加围栏：EN (4-methoxyphenyl)sulfonyl、ZH (4-甲氧基苯基)磺酰基；非 S 桥或前端为简单取代基（propan-2-yl）返回 None 走平铺融合。"""
    bridge = MONONUCLEAR_BRIDGE.get((stem_en, stem_zh))  # P-67.1.4.1.3 复合前缀：P 酰基经 O/N/S 桥连母体，整体加方括号后接桥后缀
    if bridge is not None and any((a.get("en") or "").endswith(s) for s in PHOSPHORYL_STEMS):
        return f"[{a['en']}]{bridge[0]}", f"[{a['zh']}]{bridge[1]}"
    if stem_en not in ("sulfinyl", "sulfonyl") or not a.get("paren"):
        return None
    if _SIMPLE_CHAIN_YL_RE.match(a.get("en") or ""):  # 直链 -yl 前端与桥融合：propan-2-ylsulfonyl（非 (propan-2-yl)sulfonyl）
        return None
    w_en = f"[{a['en']}]" if "(" in a["en"] else f"({a['en']})"  # 前端自带括号时升级方括号（P-16.5.2 嵌套标记）
    w_zh = f"[{a['zh']}]" if "(" in a["zh"] else f"({a['zh']})"
    return f"{w_en}{stem_en}", f"{w_zh}{stem_zh}基"  # ZH 前端基不可省（(4-甲氧基苯基)磺酰基，非 …苯磺酰基）


def _alpha_key(name: str) -> str:
    """P-14.5 字母序比较键：以小写字母为准，忽略位次/括号/连字符等非字母字符。"""
    return "".join(c for c in name.lower() if c.isalpha())


# free 母体名 → P-29 -yl 取代基形式的双语转换（原 tools/free_to_yl.py 并入）。
_MULT_OL_SUF = re.compile(r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)ol\b")  # 官能团后缀 → 前缀转换

_MONONUCLEAR_NAMES = ("oxidane", "azane", "sulfane", "sulfinyl", "sulfonyl", "imine")  # 本模块转换的单核母体氢化物（P-15.4.1 表 2.1 → 表 1.5 'a' 前缀体系，数据见 constants.MONONUCLEAR_HYDRIDES）；P 酰基（phosphoryl/phosphanyl）走 L5 的 P 专用管线，不经此表。


def _anilino_en(base: str, yl: str) -> str:
    """N-苯基（可带环取代基）的 azane 去氢：phenyl→anilino、4-chlorophenyl→4-chloroanilino（P-62.2.1.1：phenylamino = anilino*）。"""
    return base[: -len("phenyl")] + "anilino" if yl == "amino" and base.endswith("phenyl") else base + yl


def _mononuclear_en(en: str) -> str | None:
    """单核氢化物 free 名 → 去氢取代基名：ethyl-oxidane → ethyloxy、phenyl-azane → anilino。"""
    for en_suf in _MONONUCLEAR_NAMES:
        if en.endswith("-" + en_suf):
            return _anilino_en(en[: -len(en_suf) - 1], MONONUCLEAR_YL[en_suf][1])
    return None


def _mononuclear_zh(zh: str) -> str | None:
    """中文组装名去氢：乙基-氧化烷 → 乙氧基、苯基-氮烷 → 苯胺基、4-氯苯基-氮烷 → 4-氯苯胺基。"""
    for en_suf in _MONONUCLEAR_NAMES:
        zh_suf, _, zy = MONONUCLEAR_YL[en_suf]
        if zh.endswith("-" + zh_suf):
            base = zh[: -len(zh_suf) - 1]
            if zy == "氨基" and base.endswith("苯基"):
                return base[: -len("苯基")] + "苯胺基"
            return zh_bridge_root(base) + zy
    return None

def _try_fg_prefix(en: str, zh: str) -> tuple[str, str] | None:
    """尝试将官能团后缀名转换为取代基前缀形式。"""
    for en_fn, zh_fn in [
        (_mononuclear_en, _mononuclear_zh),
    ]:
        en_out = en_fn(en)
        zh_out = zh_fn(zh) if en_out else None
        if en_out and zh_out:
            return en_out, zh_out
    return None

def free_to_yl(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool]:
    """在键合位次处把 free 母体名转 P-29 -yl 双语形式（P-63.2.2 醇→烷氧基、P-63.2.1 硫醇→烷基硫基、P-62.2 1° 胺→烷基氨基）。"""
    fg = _try_fg_prefix(en, zh)
    if fg is not None:
        need_paren = fg[0].endswith("amino") and fg[0] != "amino"  # P-29.3.6：复合前缀（methylamino=CH3-NH-，非普通 amino）需括号与两个独立取代基区分。
        need_paren = need_paren or (fg[0].endswith("anilino") and fg[0] != "anilino")  # 带环取代基的 anilino（4-chloroanilino）与 …phenylamino 同理需括号（gold：(4-chloroanilino)benzoic acid）；裸 anilino 免括。
        return fg[0], fg[1], need_paren
    return None


def _mononuclear_radical_names(numbered: dict) -> tuple[str, str] | None:
    """杂原子锚点自由基：单核氢化物母体+烷基取代基经 free_to_yl 转标准名（*OCC→ethoxy；azane 双烷基按 P-62.2 字母序、同烷基 di-）；零/多取代基或名缺失返回 None 明确失败。"""
    parent = numbered.get("parent") or {}
    stem_en, stem_zh = parent.get("stem_en"), parent.get("stem_zh")
    if not stem_en or not stem_zh:
        return None
    subs = [s for s in (numbered.get("substituents") or []) if s.get("en") and s.get("zh")]
    if not subs:
        return MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    if stem_en in PHOSPHORYL_STEMS:  # P-67.1.4.1.1.5：P 酰基前缀按取代基拼接，取代基可 1–3 个（hydroxy(methoxy)phosphoryl）
        return _phosphoryl_sub_names(subs, stem_en, stem_zh, parent.get("mol"))
    if len(subs) == 1:
        a = subs[0]
        if stem_en == "sulfonyl" and a["en"].endswith("amino") and a["zh"].endswith("氨基"):  # P-66.1.1.4.2 + Glossary：(phenylamino)sulfonyl = phenylsulfamoyl*；N-取代基与 sulfamoyl 融合（丁基(甲基)sulfamoyl），不走 amino 围栏
            return a["en"][: -len("amino")] + "sulfamoyl", a["zh"] + "磺酰基"
        bridge = _bridge_enclosed_names(a, stem_en, stem_zh)  # 围栏在 L3 一次定形，中英文同步产出，L5 前缀渲染不再二次拆分
        if bridge is not None:
            if (stem_en, stem_zh) != ("azane", "氮烷"):  # N 桥复合前缀（…phosphoryl]amino）作取代基时仍须 L5 整体围栏（P-16.5.2 嵌套），O/S 桥名下自带围栏不再加
                numbered["bridge_self_enclosed"] = True
            return bridge
        if stem_en == "azane":
            amido = AMIDO_RETAINED.get(a.get("en") or "")  # P-66.1.1.4.3 方法 1：单 N-酰基（乙酰/甲酰/苯甲酰）残基收成 amido 保留式（acetamido…），不走 free_to_yl 的 acylamino 系统式；其余 R 保持方法 2。
            if amido is not None:
                return amido
            if _azane_acyl_stereo_lead(a.get("en") or "") or _azane_sub_needs_paren(a):  # 方法 2 需把内层组整体括起再加 amino（P-29.3.2 复合前缀括号，见 AZANE_PAREN_SUF）：带立体描述符的复杂酰基残基（肽类 N-酰基氨基酸）否则会与 N-端 amino 位次歧义；芳香酰基/acetyl 系见白名单；无立体简单酰与开链酰保持融合平铺。
                w_en = f"[{a['en']}]" if "(" in a["en"] else f"({a['en']})"  # 内层已含括号（立体描述符）时升级为方括号（P-16.5.2 嵌套）
                w_zh = f"[{a['zh']}]" if "(" in a["zh"] else f"({a['zh']})"
                return f"{w_en}amino", f"{w_zh}氨基"
        en, zh = free_to_yl(f"{a['en']}-{stem_en}", f"{a['zh']}-{stem_zh}", 1,
                            paren=bool(a.get("paren")))[:2]
        return _retained_alkoxy(en, zh) if stem_en == "oxidane" else (en, zh)
    zero = MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    if (stem_en, stem_zh) != ("azane", "氮烷") or zero is None or len(subs) != 2:  # 多取代基仅 N（azane）双烷基成立：O/S 双烷基非标准自由基，明确失败。
        return None
    ordered = sorted(subs, key=lambda s: s["en"])
    if len({s["en"] for s in ordered}) == 1:
        base = ordered[0]
        return (f"{MULT_EN[len(ordered)]}{base['en']}{zero[0]}",
                f"{MULT_ZH[len(ordered)]}{zh_bridge_root(base['zh'])}{zero[1]}")
    aryl = [s for s in ordered if s["en"].endswith("phenyl") and s["zh"].endswith("苯基")]
    if len(aryl) == 1:  # P-62.2.1.1：N-芳基-N-某基胺取 anilino，非芳基 N-取代基以 N- 前缀（4-fluoro-N-propan-2-ylanilino）；两前缀按 P-14.5 字母序
        ring = aryl[0]
        other = next(s for s in ordered if s is not ring)
        ring_en, ring_zh = ring["en"][: -len("phenyl")], ring["zh"][: -len("苯基")]
        if _alpha_key(ring_en) <= _alpha_key(other["en"]):
            return (f"{ring_en}-N-{other['en']}anilino" if ring_en else f"N-{other['en']}anilino",
                    f"{ring_zh}-N-{other['zh']}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
        return (f"N-{other['en']}-{ring_en}anilino" if ring_en else f"N-{other['en']}anilino",
                f"N-{other['zh']}-{ring_zh}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
    first, rest = ordered[0], ordered[1:]  # 双不同 N-取代基：字母序首基平铺，其后各基分别加括号紧贴 amino（P-62.2.2.1：多取代氨基须逐基消歧，2-chloroethylethylamino → 2-chloroethyl(ethyl)amino）。
    return (first["en"] + "".join(f"({s['en']})" for s in rest) + zero[0],
            first["zh"] + "".join(f"({s['zh']})" for s in rest) + zero[1])


def _indicated_h_prefix(parent: dict) -> str:
    """把 L4 定好的指示氢位次拼成前缀（'1H-' / '1H,2H-'），无位次则为空串（P-58.2.1）。"""
    locants = parent.get("indicated_h_locants") or ()
    return f"{','.join(f'{l}H' for l in locants)}-" if locants else ""


def _ensure_fused_stem(numbered: dict) -> bool:
    """未注册稠环词干注入：parent 无词干但有 fused_tree 时用 fused_parent_names 组装稠合 base 名（已注册词干由 L2 注入不进入）；返回 False 表示组装失败（显式 unsupported）。"""
    parent = numbered.get("parent") or {}
    if parent.get("stem_en") and parent.get("stem_zh"):
        return True
    node = parent.get("fused_tree")
    if node is None:
        return True
    mol = parent.get("mol")
    if mol is None:
        return False
    from namepredict.layer5.fused_namer import fused_parent_names
    name = fused_parent_names(mol, node)
    if name is None or not name[0] or not name[1]:
        return False
    pre = _indicated_h_prefix(parent)  # 稠合 base 名此前不带指示氢，L4 已按整体编号定好位次（P-58.2.1），统一在此补到最前端。
    parent["stem_en"], parent["stem_zh"] = pre + name[0], pre + name[1]
    return True


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """链引擎按表 kind 派发，再转具体 worker。"""
    if kind == "phosphate":  # 无机功能母体整名（P 中心无碳词干）：由 numbered 母体计数 + o_side 臂组装。
        from namepredict.layer5.phosphate import phosphate_names

        return phosphate_names(numbered)
    if kind in EXO_RING_SUF:  # 环外主基走 exocyclic worker（苯单取代→None 回落 chain_engine 保留名 variant）；开链同类无关（relation in_skeleton）。
        exo = _exocyclic_ring_names(n, numbered)
        if exo:
            return exo
    if kind == "radical" and (numbered.get("parent") or {}).get("radical_anchor_element"):
        return _mononuclear_radical_names(numbered)
    entry = _KIND_TABLE.get(kind)
    if kind == "acyl_halide":  # 酰卤后缀随实际卤素（F/Cl/Br/I）变：由 L1 检测、L2 保留的 parent.hal_z 选择 spec。
        entry = _ACYL_HALIDE_BY_HAL.get((numbered.get("parent") or {}).get("hal_z")) or entry
    if entry is not None:
        sid = _scaffold_id(numbered)
        if kind == "alkane" and (numbered.get("parent") or {}).get("fused_tree") and sid != "benzene":  # 未注册稠环无 FG: 词干注入(_ensure_fused_stem)已完成, 返回稠合 base 名(不走 chain_engine 拼 ane)。
            return _parent_stem_names(numbered)
        if sid == "benzene" and kind == "alkane":  # 苯 base：无主 FG 的苯，母体名由 sid 驱动（L2 已把纯苯 kind 收敛为 alkane）。
            return ("benzene", "苯")
        ring_stem = _ring_stem(numbered)
        if ring_stem:  # 环式 FG 的 locant omit 由 L4 的 omit 标志决定，aromatic 仅对苯环置真（苯醇→酚，杂环醇→醇）；coda 重置为空（杂环词干 pyridin/furan 已完整，不再接饱和链 "an"）。
            stem_en, stem_zh = ring_stem
            entry = replace(entry, stem=(stem_en, stem_zh), coda="",
                            omit_rule=lambda n, loc, omit: bool(omit), aromatic=(sid == "benzene"))
        elif sid == "carbocycle":  # 单环饱和烃自由基按 P-29.2 方法 1 省略自由价 1 位（cyclopentyl/cyclohexyl）；其余 FG 沿用 L4 omit 标志（cyclohexanol 等）。环烯走 unsat 段不受影响。
            rule = (lambda n, loc, omit: loc == 1) if kind == "radical" \
                else (lambda n, loc, omit: bool(omit))
            entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=rule)
        sc_variant = (entry.variant or {}).get(sid)  # 苯环单 FG → scaffold 专属保留名 variant (phenol/benzoic…); 开链取 None 键 (acid 草酸)。
        if sc_variant is not None:
            entry = replace(entry, variant=sc_variant)
        return _chain_names(entry, n, numbered)

    stem = _parent_stem_names(numbered)
    return stem

def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    """取母体已算出的双语词干（stem_en/stem_zh）。"""
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None

def _parent_n(numbered: dict) -> tuple[str | None, int]:
    """取母体 kind 与碳数 n（无则 0）。"""
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)


def _unsupported(n: int, kind: str | None) -> NameResult:
    """构造 unsupported 失败结果并携带碳数与 kind。"""
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})

# 名称拼接：前缀与母体组合（原 benzene_names.py 并入；P-22.1.3）。
def _needs_join_hyphen(body: str) -> bool:
    """词干以数字（1,3-thiazole）、方括号位次集（[1,2,4]triazolo[1,5-a]pyridine）或 1H-（1H-pyrrole）开头时，与前缀需连字符分隔。"""
    return bool(body) and (body[0].isdigit() or body[0] == "[" or body.startswith("1H-"))


def join_parent_name(prefix: str, parent: str) -> str:
    """拼接前缀与母体名（数字/1H- 前导时加连字符）。"""
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if _needs_join_hyphen(stem) else f"{prefix}{stem}"
    return f"{stereo}{body}"

_ZH_PLAIN_YL_RE = re.compile(r"^[^()\[\]\d]*-\d+-基$")  # 仅带位次的基（噻吩-3-基、丙-2-基）：中文酯名保留「基」但不加围栏（乙酸噻吩-3-基酯）。
_ZH_SIMPLE_YL_RE = re.compile(r"^[^\-()\[\]\d]+基$")  # 无位次无取代的简单烃基（乙基/苄基/叔丁基）：中文酯名习用省「基」（乙酸乙酯、十八酸苄酯）。


def _zh_alkoxy_part(name: str) -> str:
    """中文酯 O 侧基名渲染（P-65.6.3） """
    if not name:
        return name
    if _ZH_SIMPLE_YL_RE.match(name):
        return name[:-1]  # 简单烃基省「基」：乙酸乙酯、十八酸苄酯
    if _ZH_PLAIN_YL_RE.match(name):
        return name  # 只带位次：乙酸噻吩-3-基酯、乙酸丙-2-基酯
    stereo, body = _stereo_lead(name)
    if not body or body.endswith("基") is False:
        return name
    if "(" in body or "[" in body or "（" in body:
        return f"{stereo}[{body}]"
    return f"{stereo}({body})"


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str] | None:
    """拼接酯名：O 侧烷基作前缀、酸侧作主体（XX 酸 YY 酯）；多酯按 P-16.3.2 用二/双(bis) 倍增（草酸二乙酯 / 己二酸双(2-乙基己基)酯）；无 O 侧取代基时输出 bare 酸酯名（benzoate/苯甲酸酯）。"""
    en, zh = names
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    alk_en, alk_zh = "", ""
    if len(o) == 1:
        alk_en, alk_zh = o[0].get("en") or "", _zh_alkoxy_part(o[0].get("zh") or "")
    elif len(o) > 1:
        names_en = [s.get("en") or "" for s in o]
        names_zh = [_zh_alkoxy_part(s.get("zh") or "") for s in o]
        if len(set(names_en)) == 1 and names_en[0]:  # 同名臂：di/tri 或 bis（复合前缀，P-16.3.2 须用 bis）
            if any(s.get("paren") for s in o):
                alk_en = f"{BIS_EN.get(len(o), '')}({names_en[0]})"
                alk_zh = f"{MULT_ZH.get(len(o), '')}{names_zh[0]}"
            else:
                m_en, m_zh = MULT_EN.get(len(o)), MULT_ZH.get(len(o))
                if not m_en or not m_zh:
                    return None
                alk_en, alk_zh = f"{m_en}{names_en[0]}", f"{m_zh}{names_zh[0]}"
        else:  # 异名臂：依次平铺（methyl ethyl oxalate）
            alk_en, alk_zh = " ".join(names_en), "".join(names_zh)
    st, body = _stereo_lead(en)
    mid = f"{pre_en}-{body}" if pre_en and _needs_join_hyphen(body) else f"{pre_en}{body}" if pre_en else body
    en = f"{alk_en} {st}{mid}" if alk_en else f"{st}{mid}"
    stz, bodyz = _stereo_lead(zh)
    midz = (f"{stz}{pre_zh}-{bodyz}" if pre_zh and _needs_join_hyphen(bodyz)
            else f"{stz}{pre_zh}{bodyz}" if pre_zh else f"{stz}{bodyz}")
    zh = f"{midz}{alk_zh}酯"
    return en, zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str] | None:
    """按 kind 分派：酯走酯拼接，其余走普通母体拼接。"""
    if kind in ("ester"):
        return join_ester_name(pre[0], pre[1], names, numbered)
    en = join_parent_name(pre[0], names[0])
    zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    return en, zh


def zh_1h_parent(en_parent: str, zh_parent: str, prefix: str) -> str:
    """环被取代时，中文保留名 1H- 母体补加 1H- 前缀。"""
    if not prefix or not en_parent.startswith("1H-") or zh_parent.startswith("1H-"):
        return zh_parent
    return f"1H-{zh_parent}"


def join_hydro_prefix(names: tuple[str, str], numbered: dict) -> tuple[str, str]:
    """把动态指示氢与 hydro 前缀依次拼到母体名前：顺序为 hydro + 指示氢 + 母体名；母体名已带静态 1H-（保留名）时不重复。"""
    parent = numbered.get("parent") or {}
    pre = parent.get("hydro_prefix") or ("", "")
    if not pre[0] and not parent.get("indicated_h_forced"):  # forced = 指示氢来自「保留母体名未隐含的芳香位 H」（P-58.2.1），无 hydro 前缀也须注入；其余动态指示氢仍只在氢化衍生物（有 hydro 前缀）时注入，否则 [nH] 互变异构型会误产 1H-pyridine。
        return names
    en, zh = names
    ind = _indicated_h_prefix(parent)
    if ind and not en.startswith(ind):
        en, zh = ind + en, ind + zh
    return (join_parent_name(pre[0], en), join_parent_name(pre[1], zh)) if pre[0] else (en, zh)


def join_ring_cation_suffix(numbered: dict, names: tuple[str, str]) -> tuple[str, str]:
    """净正电荷分子的环内 N+/O+ → 母体名缀 `-{位次}-ium`（P-62.4.1；chromene → chromenylium）。非该情形原样返回。"""
    parent = numbered.get("parent") or {}
    mol = parent.get("mol")
    chain = parent.get("chain") or []
    if mol is None or not chain:
        return names
    if sum(a.GetFormalCharge() for a in mol.GetAtoms()) < 0:
        return names  # 净正/中（盐、两性离子）都可能含环阳离子；净负分子不处理
    en, zh = names
    if "ium" in en:
        return names
    charged = [a.GetIdx() for a in mol.GetAtoms()
               if a.GetFormalCharge() == 1 and a.GetSymbol() in ("N", "O")
               and a.IsInRing() and a.GetIdx() in chain]
    if not charged:
        return names
    stem_en = parent.get("stem_en") or ""
    if not stem_en:
        return names
    labels = (parent.get("numbering_scaffold") or {}).get("labels")
    idx = chain.index(charged[0])
    loc = labels[idx] if labels and len(labels) == len(chain) else str(idx + 1)
    if stem_en.endswith("ene"):  # 色烯型氧鎓保留名：chromene → chromenylium（位次隐含；取代基中词干已省 e 为 chromen-…）
        base, ium = stem_en[:-1], f"{stem_en[:-3]}enylium"
    elif stem_en.endswith("e"):
        base, ium = stem_en[:-1], f"{stem_en[:-1]}-{loc}-ium"
    else:
        base, ium = stem_en, f"{stem_en}-{loc}-ium"
    if stem_en in en:  # 完整母体名：整词干替换
        return en.replace(stem_en, ium, 1), zh
    token = base + "e" if base + "e" in en else base
    if token in en:
        if stem_en.endswith("ene"):  # 色烯型氧鎓：chromene/chromen → chromenylium
            return en.replace(token, ium, 1), zh
        at = en.index(token) + len(token)  # 其余在词干后插入 -{位次}-ium
        return en[:at] + f"-{loc}-ium" + en[at:], zh
    return names

def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    """组装入口：取名 → 前缀 → 阴离子/R-S/金属盐后缀。"""
    from namepredict.layer5.stereo import join_rs_prefix
    kind, n = _parent_n(numbered)
    if not _ensure_fused_stem(numbered):
        return _unsupported(n, kind)
    names = _names_for(kind, n, numbered)#n 碳数
    if not names:
        return _unsupported(n, kind)
    names = join_hydro_prefix(names, numbered)
    names = join_ring_cation_suffix(numbered, names)

    joined = join_kind_name(kind, _prefix_for(numbered, kind, n), names, numbered)
    if joined is None:
        return _unsupported(n, kind)
    en, zh = joined
    en, zh = join_anion_names(numbered, en, zh)
    en, zh = join_rs_prefix(numbered, en, zh)
    en, zh = join_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
