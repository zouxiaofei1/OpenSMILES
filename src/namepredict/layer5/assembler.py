"""L5 名称组装主入口：取名后拼前缀、立体与盐后缀。"""
from __future__ import annotations
from dataclasses import replace
import re
from namepredict.constants import (
    ALKOXY_YLOXY_EN, ALKOXY_YLOXY_ZH, AMIDO_RETAINED, AZANE_PAREN_SUF, BIS_EN, BRIDGE_YL_SUFFIX,
    BRIDGE_ZH_YL_SUFFIX, ESTER_O_SIDE_KINDS,
    MONONUCLEAR_BRIDGE, MONONUCLEAR_YL, MONONUCLEAR_ZERO_YL, MULT_EN,
    MULT_ZH, PHOSPHORYL_STEMS, ZH_DIGITS, zh_bridge_root,
)
from namepredict.layer5.chain_engine import (
    _ACYL_HALIDE_BY_HAL, _BENZENE_RETAINED, _KIND_TABLE, _chain_names,
)
from namepredict.layer5.stems import (
    _metal_en_prefix, _metal_zh_suffix, join_anion_names,
)
from namepredict.layer5.assembler_prefixes import _SIMPLE_CHAIN_YL_RE, _enclose, _prefix_for
from namepredict.layer5.stereo import _split_stereo_lead as _stereo_lead
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult:
    """构造失败 NameResult（success=False，meta 供诊断）。"""
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    """构造成功 NameResult（success=True，记录耗时与来源）。"""
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

def _bracket_bridge_suffix(en: str, zh: str) -> tuple[str, str]:
    """复合组分加方括号；-yl 型烷氧/硫基把桥后缀挪到括号外（P-16.5.2）。"""
    out: list[str] = []
    for text, table in ((en, BRIDGE_YL_SUFFIX), (zh, BRIDGE_ZH_YL_SUFFIX)):
        for yl_suf, bridge in table:
            if text.endswith(yl_suf):
                text = f"[{text[: -len(bridge)]}]{bridge}"
                break
        else:
            text = f"[{text}]"
        out.append(text)
    return out[0], out[1]


def _oxido_arm(s: dict, mol) -> tuple[str, str] | None:
    """P 上氧负离子臂的取代基名 oxido/氧化（P-72.6.2）。"""
    atoms = s.get("atoms") or []
    if mol is None or len(atoms) != 1:
        return None
    a = mol.GetAtomWithIdx(int(atoms[0]))
    if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1:
        return ("oxido", "氧化")
    return None


def _phosphoryl_sub_names(subs: list[dict], stem_en: str, stem_zh: str, mol=None) -> tuple[str, str] | None:
    """P 酰基前缀的取代基拼接（P-67.1.4.1.1.5，简单基平铺/括起）。"""
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
    if len(rows) == 1 and rows[0][2] > 1:  # 同基倍增（P-16.3.2 简单基用 di-）
        en, zh, m = rows[0]
        m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
        if not m_en or not m_zh:
            return None
        return f"{m_en}{en}{stem_en}", f"{m_zh}{zh}{stem_zh}基"
    compound = any("(" in en or "[" in en for en, _, _ in rows)  # 含自身带括号的复合组分：逐组分连字符 + 方括号围栏（P-16.5.2）
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
    """单 N-酰基残基名是否带立体前导且为酰基词干。"""
    if not en.startswith("("):
        return False
    tag, stem = _stereo_lead(en)
    return bool(tag) and (stem.endswith("oyl") or "carbonyl" in stem)

def _azane_sub_needs_paren(a: dict) -> bool:
    """azane 单取代基是否需括起再缀 amino（P-16.5.1.1）。"""
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
    """S 桥前端为复合取代基时加围栏，否则返回 None。"""
    bridge = MONONUCLEAR_BRIDGE.get((stem_en, stem_zh))  # P-67.1.4.1.3 P 酰基经 O/N/S 桥连母体加方括号
    if bridge is not None and any((a.get("en") or "").endswith(s) for s in PHOSPHORYL_STEMS):
        return f"[{a['en']}]{bridge[0]}", f"[{a['zh']}]{bridge[1]}"
    if stem_en not in ("sulfinyl", "sulfonyl") or not a.get("paren"):
        return None
    if _SIMPLE_CHAIN_YL_RE.match(a.get("en") or ""):  # 直链 -yl 前端与桥融合（propan-2-ylsulfonyl）
        return None
    w_en = _enclose(a["en"])  # 前端自带括号时升级方括号（P-16.5.2 嵌套标记）
    w_zh = _enclose(a["zh"])
    return f"{w_en}{stem_en}", f"{w_zh}{stem_zh}基"  # ZH 前端基不可省（(4-甲氧基苯基)磺酰基，非 …苯磺酰基）


def _alpha_key(name: str) -> str:
    """P-14.5 字母序比较键：以小写字母为准，忽略位次/括号/连字符等非字母字符。"""
    return "".join(c for c in name.lower() if c.isalpha())


# free 母体名 → P-29 -yl 取代基形式的双语转换。
_MONONUCLEAR_NAMES = ("oxidane", "azane", "sulfane", "sulfinyl", "sulfonyl", "imine")  # 本模块转换的单核母体氢化物（P-15.4.1 表 2.1）。


def _anilino_en(base: str, yl: str) -> str:
    """N-苯基（可带环取代基）的 azane 去氢为 anilino。"""
    return base[: -len("phenyl")] + "anilino" if yl == "amino" and base.endswith("phenyl") else base + yl


def _mononuclear_en(en: str) -> str | None:
    """单核氢化物 free 名 → 去氢取代基名。"""
    for en_suf in _MONONUCLEAR_NAMES:
        if en.endswith("-" + en_suf):
            return _anilino_en(en[: -len(en_suf) - 1], MONONUCLEAR_YL[en_suf][1])
    return None


def _mononuclear_zh(zh: str) -> str | None:
    """中文组装名去氢：乙基-氧化烷 → 乙氧基等。"""
    for en_suf in _MONONUCLEAR_NAMES:
        zh_suf, _, zy = MONONUCLEAR_YL[en_suf]
        if zh.endswith("-" + zh_suf):
            base = zh[: -len(zh_suf) - 1]
            if zy == "氨基" and base.endswith("苯基"):
                return base[: -len("苯基")] + "苯胺基"
            return zh_bridge_root(base) + zy
    return None

def _fg_prefix(en: str, zh: str) -> tuple[str, str] | None:
    """尝试将单核氢化物母体名转换为取代基前缀形式（双语须同时命中）。"""
    en_out = _mononuclear_en(en)
    zh_out = _mononuclear_zh(zh) if en_out else None
    return (en_out, zh_out) if en_out and zh_out else None

def free_to_yl(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool] | None:
    """把 free 母体名转 P-29 -yl 双语形式（P-63.2.2 醇/胺）。"""
    fg = _fg_prefix(en, zh)
    if fg is not None:
        need_paren = fg[0].endswith("amino") and fg[0] != "amino"  # P-29.3.6：复合前缀（methylamino）需括号区分取代基
        need_paren = need_paren or (fg[0].endswith("anilino") and fg[0] != "anilino")  # 带环取代基的 anilino 同理需括号，裸 anilino 免括。
        return fg[0], fg[1], need_paren
    return None


def _mononuclear_radical_names(numbered: dict) -> tuple[str, str] | None:
    """杂原子锚点自由基：经 free_to_yl 转标准名（P-62.2）。"""
    parent = numbered.get("parent") or {}
    stem_en, stem_zh = parent.get("stem_en"), parent.get("stem_zh")
    if not stem_en or not stem_zh:
        return None
    subs = [s for s in (numbered.get("substituents") or []) if s.get("en") and s.get("zh")]
    if not subs:
        return MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    if stem_en in PHOSPHORYL_STEMS:  # P-67.1.4.1.1.5：P 酰基前缀按取代基拼接
        return _phosphoryl_sub_names(subs, stem_en, stem_zh, parent.get("mol"))
    if len(subs) == 1:
        a = subs[0]
        if stem_en == "sulfonyl" and a["en"].endswith("amino") and a["zh"].endswith("氨基"):  # P-66.1.1.4.2：N-取代基与 sulfamoyl 融合
            return a["en"][: -len("amino")] + "sulfamoyl", a["zh"] + "磺酰基"
        bridge = _bridge_enclosed_names(a, stem_en, stem_zh)  # 围栏在 L3 一次定形，中英文同步产出，L5 前缀渲染不再二次拆分
        if bridge is not None:
            if (stem_en, stem_zh) != ("azane", "氮烷"):  # N 桥复合前缀仍须 L5 整体围栏（P-16.5.2）
                numbered["bridge_self_enclosed"] = True
            return bridge
        if stem_en == "azane":
            amido = AMIDO_RETAINED.get(a.get("en") or "")  # P-66.1.1.4.3 方法 1：单 N-酰基残基收成 amido 保留式
            if amido is not None:
                return amido
            if _azane_acyl_stereo_lead(a.get("en") or "") or _azane_sub_needs_paren(a):  # 方法 2 需把内层组整体括起再加 amino（P-29.3.2）
                w_en = _enclose(a["en"])  # 内层已含括号（立体描述符）时升级为方括号（P-16.5.2 嵌套）
                w_zh = _enclose(a["zh"])
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
    if len(aryl) == 1:  # P-62.2.1.1：N-芳基-N-某基胺取 anilino
        ring = aryl[0]
        other = next(s for s in ordered if s is not ring)
        ring_en, ring_zh = ring["en"][: -len("phenyl")], ring["zh"][: -len("苯基")]
        if _alpha_key(ring_en) <= _alpha_key(other["en"]):
            return (f"{ring_en}-N-{other['en']}anilino" if ring_en else f"N-{other['en']}anilino",
                    f"{ring_zh}-N-{other['zh']}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
        return (f"N-{other['en']}-{ring_en}anilino" if ring_en else f"N-{other['en']}anilino",
                f"N-{other['zh']}-{ring_zh}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
    first, rest = ordered[0], ordered[1:]  # 双不同 N-取代基：首基平铺，其余各基加括号（P-62.2.2.1）
    return (first["en"] + "".join(f"({s['en']})" for s in rest) + zero[0],
            first["zh"] + "".join(f"({s['zh']})" for s in rest) + zero[1])


def _indicated_h_prefix(parent: dict) -> str:
    """把 L4 定好的指示氢位次拼成前缀（P-58.2.1）。"""
    locants = parent.get("indicated_h_locants") or ()
    return f"{','.join(f'{l}H' for l in locants)}-" if locants else ""


def _ensure_fused_stem(numbered: dict) -> bool:
    """未注册稠环词干注入：由 fused_tree 组装稠合 base 名。"""
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
    pre = _indicated_h_prefix(parent)  # L4 已按整体编号定好指示氢位次（P-58.2.1），此处补到最前
    parent["stem_en"], parent["stem_zh"] = pre + name[0], pre + name[1]
    return True

def _phosphate_arm_zh(zh: str) -> str:
    """磷酸臂中文词：简单基去「基」，复合名原样保留。"""
    if not zh.endswith("基") or "-" in zh or zh.startswith("("):
        return zh
    stem = zh[:-1]
    if len(stem) >= 2 and all(c in ZH_DIGITS + "十" for c in stem):
        return f"{stem}烷基"
    return stem


def join_phosphate_name(names: tuple[str, str], numbered: dict) -> tuple[str, str] | None:
    """磷酸整名（P-67.1.3）：O-侧臂 + 词尾 + 金属盐。"""
    tail_en, tail_zh = names
    arms = _join_o_side_arms(_o_side_arms(numbered), group=True, arm_zh_fn=_phosphate_arm_zh)
    if arms is None:
        return None
    alk_en, alk_zh = arms
    parent = numbered.get("parent") or {}
    salt_meta = parent.get("salt_meta") or {}
    metal_en = _metal_en_prefix(salt_meta)
    metal_zh = _metal_zh_suffix(salt_meta)
    if not alk_en:  # 无 O-侧臂
        if metal_en and metal_zh:  # 酸式盐：金属直接缀在酸式词后（磷酸二氢钾），不带"酯"
            return f"{metal_en} {tail_en}", f"{tail_zh}{metal_zh}"
        if int(parent.get("n_om") or 0) > 0:  # 游离磷酸根：负电荷不标注，用基本根词
            return tail_en, f"{tail_zh}根"
        return tail_en, tail_zh
    if metal_en and metal_zh:  # 酯盐：金属名前置，中文"{金属}盐"后置
        return f"{metal_en} {alk_en} {tail_en}", f"{tail_zh}{alk_zh}酯 {metal_zh}盐"
    return f"{alk_en} {tail_en}", f"{tail_zh}{alk_zh}酯"


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """链引擎按表 kind 派发，再转具体 worker。"""
    parent = numbered.get("parent") or {}
    if kind == "radical" and parent.get("radical_anchor_element"):
        return _mononuclear_radical_names(numbered)
    entry = _KIND_TABLE.get(kind)
    if kind == "acyl_halide":
        entry = _ACYL_HALIDE_BY_HAL.get(parent.get("hal_z")) or entry
    if entry is not None:
        sid = parent.get("scaffold_id")
        if kind == "alkane" and parent.get("fused_tree") and sid != "benzene":  # 未注册稠环无 FG：词干注入已完成，返回稠合 base 名
            return _parent_stem_names(numbered)
        if sid == "benzene" and kind == "alkane":  # 苯 base：母体名由 sid 驱动
            return ("benzene", "苯")
        if parent.get("stem_en") and parent.get("stem_zh"):  # 环式 FG 的 locant omit 由 L4 决定，coda 置空
            entry = replace(entry, stem=(parent["stem_en"], parent["stem_zh"]), coda="",
                            omit_rule=lambda n, loc, omit: bool(omit), aromatic=(sid == "benzene"))
        elif sid == "carbocycle":  # 单环饱和烃自由基按 P-29.2 省略 1 位自由价，余沿用 L4 omit
            rule = (lambda n, loc, omit: loc == 1) if kind == "radical" \
                else (lambda n, loc, omit: bool(omit))
            entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=rule)
        # 苯环单 FG 取 scaffold 专属保留名（P-61.2 表）。
        sc_variant = _BENZENE_RETAINED[kind] if sid == "benzene" and kind in _BENZENE_RETAINED \
            else (entry.variant or {}).get(sid)
        if sc_variant is not None:
            entry = replace(entry, variant=sc_variant)
        result = _chain_names(entry, n, numbered)
        print("\n\nresult:   ",result)
        return result

    return _parent_stem_names(numbered)

def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    """取母体双语词干；结尾 'e' 省略按 P-60.2(a) 判定。"""
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None

def _unsupported(n: int, kind: str | None) -> NameResult:
    """构造 unsupported 失败结果并携带碳数与 kind。"""
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})

# 名称拼接：前缀与母体组合（P-22.1.3）。
def _needs_join_hyphen(body: str) -> bool:
    """词干以数字、方括号位次集或 1H- 开头时需与前缀连字符分隔。"""
    return bool(body) and (body[0].isdigit() or body[0] == "[" or body.startswith("1H-"))


def join_parent_name(prefix: str, parent: str) -> str:
    """拼接前缀与母体名（数字/1H- 前导时加连字符）。"""
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if _needs_join_hyphen(stem) else f"{prefix}{stem}"
    return f"{stereo}{body}"

_ZH_PLAIN_YL_RE = re.compile(r"^[^()\[\]\d]*-\d+-基$")  # 仅带位次的基：中文酯名保留「基」不加围栏
_ZH_SIMPLE_YL_RE = re.compile(r"^[^\-()\[\]\d]+基$")  # 简单烃基：中文酯名习用省「基」


def _zh_alkoxy_part(name: str) -> str:
    """中文酯 O 侧基名渲染（P-65.6.3） """
    if not name:
        return name
    if _ZH_SIMPLE_YL_RE.match(name):
        return name[:-1]  # 简单烃基省「基」：乙酸乙酯、十八酸苄酯
    if _ZH_PLAIN_YL_RE.match(name):
        return name  # 只带位次：乙酸噻吩-3-基酯、乙酸丙-2-基酯
    stereo, body = _stereo_lead(name)
    if not body or not body.endswith("基"):
        return name
    # 前端已含围栏（圆/方/全角）时升级方括号，否则加圆括号。
    enclose = "(" in body or "[" in body or "（" in body
    return f"{stereo}[{body}]" if enclose else f"{stereo}({body})"


def _o_side_arms(numbered: dict) -> list[dict]:
    """母体的 O-侧臂（酯/磷酸的烷氧基臂）取代基列表，非 O-侧母体为空。"""
    return [s for s in (numbered.get("substituents") or []) if s.get("o_side")]


def _join_o_side_arms(arms: list[dict], *, group: bool, arm_zh_fn) -> tuple[str, str] | None:
    """O-侧臂双语拼接"""
    if not arms:
        return "", ""
    if len(arms) == 1:
        return arms[0].get("en") or "", arm_zh_fn(arms[0].get("zh") or "")
    if group:
        table: dict[str, list] = {}
        for s in arms:
            en = (s.get("en") or "").strip()
            if not en:
                continue
            row = table.setdefault(en, [en, arm_zh_fn(s.get("zh") or ""), 0])
            row[2] += 1
        parts_en: list[str] = []
        parts_zh: list[str] = []
        for en in sorted(table):
            _, zh, m = table[en]
            if m == 1:
                parts_en.append(en)
                parts_zh.append(zh)
                continue
            m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
            if not m_en or not m_zh:
                return None
            parts_en.append(f"{m_en}{en}")
            parts_zh.append(f"{m_zh}{zh}")
        return " ".join(parts_en), "".join(parts_zh)
    names_en = [s.get("en") or "" for s in arms]
    names_zh = [arm_zh_fn(s.get("zh") or "") for s in arms]
    if len(set(names_en)) == 1 and names_en[0]:  # 同名臂：di/tri 或 bis（复合前缀，P-16.3.2 须用 bis）
        if any(s.get("paren") for s in arms):
            return (f"{BIS_EN.get(len(arms), '')}({names_en[0]})",
                    f"{MULT_ZH.get(len(arms), '')}{names_zh[0]}")
        m_en, m_zh = MULT_EN.get(len(arms)), MULT_ZH.get(len(arms))
        if not m_en or not m_zh:
            return None
        return f"{m_en}{names_en[0]}", f"{m_zh}{names_zh[0]}"
    return " ".join(names_en), "".join(names_zh)  # 异名臂：依次平铺（methyl ethyl oxalate）


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str] | None:
    """拼接酯名：O 侧作前缀、酸侧作主体（P-16.3.2 倍增）。"""
    en, zh = names
    arms = _join_o_side_arms(_o_side_arms(numbered), group=False, arm_zh_fn=_zh_alkoxy_part)
    # print("_join_o_side_arms",arms,numbered)
    if arms is None:
        return None
    alk_en, alk_zh = arms
    mid = join_parent_name(pre_en, en)
    en = f"{alk_en} {mid}" if alk_en else mid
    zh = f"{join_parent_name(pre_zh, zh)}{alk_zh}酯"
    return en, zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str] | None:
    """按 kind 分派：O-侧臂母体（酯/磷酸）走 O-侧拼接，其余走普通母体拼接。"""
    if kind in ESTER_O_SIDE_KINDS:
        if kind == "phosphate":
            return join_phosphate_name(names, numbered)
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
    """把 hydro 与指示氢前缀依次拼到母体名前。"""
    parent = numbered.get("parent") or {}
    pre = parent.get("hydro_prefix") or ("", "")
    if not pre[0] and not parent.get("indicated_h_forced"):  # forced = 指示氢来自保留母体名未隐含的芳香位 H（P-58.2.1）
        return names
    en, zh = names
    ind = _indicated_h_prefix(parent)
    if ind and not en.startswith(ind):
        en, zh = ind + en, ind + zh
    return (join_parent_name(pre[0], en), join_parent_name(pre[1], zh)) if pre[0] else (en, zh)


def join_ring_cation_suffix(numbered: dict, names: tuple[str, str]) -> tuple[str, str]:
    """环内 N+/O+ → 母体名缀 -{位次}-ium（P-62.4.1）。"""
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
    if stem_en.endswith("ene"):  # 色烯型氧鎓保留名：chromene → chromenylium
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
    from namepredict.layer5.stereo import join_ez_prefix, join_rs_prefix
    parent = numbered.get("parent") or {}  # 母体 kind 与碳数 n（无则 0）
    kind, n = parent.get("kind"), int(parent.get("n_carbons") or 0)
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
    en, zh = join_ez_prefix(numbered, en, zh)   # 母体外挂双键的 E/Z 先补，再由 R/S 归并排序
    en, zh = join_rs_prefix(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
