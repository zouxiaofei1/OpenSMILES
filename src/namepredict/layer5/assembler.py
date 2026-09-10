"""L5 名称组装主入口：由链引擎取名后按 kind 拼接前缀、立体（E/Z、R/S）与盐类后缀。"""
from __future__ import annotations
from dataclasses import replace

from namepredict.constants import AMIDO_RETAINED, MULT_EN, MULT_ZH
from namepredict.layer5.chain_engine import _ACYL_HALIDE_BY_HAL, _KIND_TABLE, _alkane_names, _chain_names
from namepredict.layer5.stems import maybe_anion_names, maybe_metal_salt_names
from namepredict.layer5.assembler_prefixes import _prefix_for
from namepredict.layer5.stereo import _split_stereo_lead as _stereo_lead
from namepredict.tools.free_to_yl import free_to_yl
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult:
    """构造失败 NameResult（success=False，meta 供诊断）。"""
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    """构造成功 NameResult（success=True，记录耗时与来源）。"""
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

# 稠环/杂环词干（en_stem, zh_stem）：由 L2 注入 parent 的 stem_en/stem_zh 派生。
# IUPAC 词干 = 母体名去尾部 e（benzene 例外，保留完整名）；aromatic 恒 True（保留母体均芳香）。
def _ring_stem(numbered: dict) -> tuple[str, str] | None:
    """保留 scaffold 的完整 IUPAC 词干（用于 -ol/-diol/-amine 等 FG 后缀拼接）；结尾 'e' 的省略交给 chain_engine._elide_parent_e 按后缀首字母判断（P-60.2(a)），避免 oxolane-3,4-diol 被错拼成 oxolan-3,4-diol。"""
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None


# 稠环/杂环完整 base 名（-carboxylic acid 用完整词干，如 naphthalene-1-carboxylic acid）。
def _ring_base(numbered: dict) -> tuple[str, str] | None:
    """保留 scaffold 的完整母体名（用于 -carboxylic acid）。"""
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
    ene = numbered.get("ene_locant")
    enes = numbered.get("ene_locants")
    if not (ene or enes):
        return f"cyclo{base[0]}", f"环{base[1]}", False
    en_core = base[0][:-3] if base[0].endswith("ane") else base[0]   # hexane → hex
    zh_core = base[1][:-1] if base[1].endswith("烷") else base[1]    # 己烷 → 己
    if enes and len(enes) >= 2:
        loc = ",".join(str(x) for x in enes)
        m_en, m_zh = MULT_EN.get(len(enes)), MULT_ZH.get(len(enes))
        if not m_en or not m_zh:
            return None, None, True
        return (f"cyclo{en_core}a-{loc}-{m_en}ene", f"环{zh_core}-{loc}-{m_zh}烯", True)
    if ene == 1:
        return f"cyclo{en_core}ene", f"环{zh_core}-1-烯", True
    return f"cyclo{en_core}-{ene}-ene", f"环{zh_core}-{ene}-烯", True


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


# 环外主基系统名后缀表（group_class → 后缀规格），六个环外 worker 共用 `_exocyclic_ring_names` 一条管线：
#   singular          (en, zh)          单取代后缀
#   plural            (en, zh) | None   多取代后缀基底（前拼 MULT_EN/MULT_ZH 倍数词）；None = 该主基无多取代系统名
#   plural_needs_loc  bool              多取代时位次缺失即放弃（P-65.2.2：多羧酸位次必带）
_EXO_SUF: dict[str, tuple] = {
    "acid":     (("carboxylic acid", "羧酸"), ("carboxylic acid", "羧酸"), True),
    "aldehyde": (("carbaldehyde", "甲醛"),    ("carbaldehyde", "甲醛"),    False),
    "ester":    (("carboxylate", "羧酸"),     None,                        False),
    "amide":    (("carboxamide", "甲酰胺"),   None,                        False),
    "nitrile":  (("carbonitrile", "甲腈"),    None,                        False),
    "acyl":     (("carbonyl", "羰基"),        None,                        False),
}


def _exocyclic_ring_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """环外主基（羧基/醛/酯/酰胺/腈/酰基头）系统名：环母体词干 + `_EXO_SUF` 后缀组装（P-65.1.7.2 酰基头 / P-65.2.2 多羧酸 / P-66.6.1.1.3 环醛）；苯单取代走 chain_engine 保留名（benzoic acid/benzaldehyde/benzonitrile…）时返回 None 回落。"""
    parent = numbered.get("parent") or {}
    facts = parent.get("principal_expression_facts")
    sid = parent.get("scaffold_id")
    spec = _EXO_SUF.get(facts.group_class.value) if facts else None
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
        # 单取代：环烯使编号不再唯一、或环上另带前缀取代（必带位次）时主基 locant 1 不可省略（P-65.2.2.1/P-66.6.1）；多取代恒带位次。
        use_loc = bool(loc) and (mult > 1 or has_unsat or _ring_extra_prefix_located(numbered))
        if use_loc:
            return (f"{en_ring}-{loc}-{suf_en}", f"{zh_ring}-{loc}-{suf_zh}")
        return (f"{en_ring}{suf_en}", f"{zh_ring}{suf_zh}")
    base = _ring_base(numbered)
    if base:  # 主基位次取 L4 已算出的 FG locant；无则省略（1 位隐含）。
        if loc:
            return (f"{base[0]}-{loc}-{suf_en}", f"{base[1]}-{loc}-{suf_zh}")
        return (f"{base[0]}{suf_en}", f"{base[1]}{suf_zh}")
    return None


_MONONUCLEAR_ZERO_YL = {
    ("oxidane", "氧化烷"): ("hydroxy", "羟基"),
    ("azane", "氮烷"): ("amino", "氨基"),
    ("sulfane", "硫烷"): ("sulfanyl", "硫基"),
    ("sulfinyl", "亚磺酰"): ("sulfinyl", "亚磺酰基"),
    ("sulfonyl", "磺酰"): ("sulfonyl", "磺酰基"),
    ("imine", "亚胺"): ("imino", "亚氨基"),
}


def _azane_acyl_stereo_lead(en: str) -> bool:
    """单 N-酰基残基名是否带前导立体描述符且为酰基词干（azane 方法 2 需括起 acyl 再缀 amino）。"""
    if not en.startswith("("):
        return False
    tag, stem = _stereo_lead(en)
    return bool(tag) and (stem.endswith("oyl") or "carbonyl" in stem)

# O 锚点自由基 -yloxy 非保留名 → IUPAC 保留烷氧基（P-66.5.2.1.2：ethoxy/propoxy/butoxy/phenoxy）。
# 尾部收拢使带取代基链也命中：2-methoxyethyloxy → 2-methoxyethoxy、3-chlorophenyloxy → 3-chlorophenoxy。
_ALKOXY_YLOXY_EN = (
    ("ethyloxy", "ethoxy"), ("propyloxy", "propoxy"), ("butyloxy", "butoxy"),
    ("phenyloxy", "phenoxy"),
)
_ALKOXY_YLOXY_ZH = (
    ("乙基氧基", "乙氧基"), ("丙基氧基", "丙氧基"), ("丁基氧基", "丁氧基"),
    ("苯基氧基", "苯氧基"),
)


def _retained_alkoxy(en: str, zh: str) -> tuple[str, str]:
    """O 锚点 -yloxy 尾部收拢为 IUPAC 保留烷氧基（未命中原样返回）。"""
    for suf_en, kept_en in _ALKOXY_YLOXY_EN:
        if en.endswith(suf_en):
            en = en[: -len(suf_en)] + kept_en
            break
    for suf_zh, kept_zh in _ALKOXY_YLOXY_ZH:
        if zh.endswith(suf_zh):
            zh = zh[: -len(suf_zh)] + kept_zh
            break
    return (en, zh)


def _mononuclear_radical_names(numbered: dict) -> tuple[str, str] | None:
    """杂原子锚点自由基：单核氢化物母体+烷基取代基经 free_to_yl 转标准名（*OCC→ethoxy；azane 双烷基按 P-62.2 字母序、同烷基 di-）；零/多取代基或名缺失返回 None 明确失败。"""
    parent = numbered.get("parent") or {}
    stem_en, stem_zh = parent.get("stem_en"), parent.get("stem_zh")
    if not stem_en or not stem_zh:
        return None
    subs = [s for s in (numbered.get("substituents") or []) if s.get("en") and s.get("zh")]
    if len(subs) == 0:
        return _MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    if len(subs) == 1:
        a = subs[0]
        if stem_en == "azane":
            amido = AMIDO_RETAINED.get(a.get("en") or "")  # P-66.1.1.4.3 方法 1：单 N-酰基（乙酰/甲酰/苯甲酰）残基收成 amido 保留式（acetamido…），不走 free_to_yl 的 acylamino 系统式；其余 R 保持方法 2。
            if amido is not None:
                return amido
            if _azane_acyl_stereo_lead(a.get("en") or ""):  # 带立体描述符的复杂酰基残基（肽类 N-酰基氨基酸）：方法 2 需把酰基名整体括起再加 amino（[…propanoyl]amino，P-29.3.2 复合前缀括号），否则融合式会与 N-端 amino 位次歧义；无立体简单酰保持融合。
                return f"[{a['en']}]amino", f"[{a['zh']}]氨基"
        en, zh = free_to_yl(f"{a['en']}-{stem_en}", f"{a['zh']}-{stem_zh}", 1,
                            paren=bool(a.get("paren")))[:2]
        return _retained_alkoxy(en, zh) if stem_en == "oxidane" else (en, zh)
    zero = _MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    if (stem_en, stem_zh) != ("azane", "氮烷") or zero is None or len(subs) != 2:  # 多取代基仅 N（azane）双烷基成立：O/S 双烷基非标准自由基，明确失败。
        return None
    ordered = sorted(subs, key=lambda s: s["en"])
    if len({s["en"] for s in ordered}) == 1:
        base = ordered[0]
        return (f"{MULT_EN[len(ordered)]}{base['en']}{zero[0]}",
                f"{MULT_ZH[len(ordered)]}{base['zh']}{zero[1]}")
    first, rest = ordered[0], ordered[1:]  # 双不同 N-取代基：字母序首基平铺，其后各基分别加括号紧贴 amino（P-62.2.2.1：多取代氨基须逐基消歧，2-chloroethylethylamino → 2-chloroethyl(ethyl)amino）。
    return (first["en"] + "".join(f"({s['en']})" for s in rest) + zero[0],
            first["zh"] + "".join(f"({s['zh']})" for s in rest) + zero[1])


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
    pre = parent.get("indicated_h") or ""  # L4 已按整体编号算好指示氢前缀（P-58.2.1），稠合名此前不带，统一在此补到最前端。
    parent["stem_en"], parent["stem_zh"] = pre + name[0], pre + name[1]
    return True


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """链引擎按表 kind 派发，再转具体 worker。"""
    if kind == "phosphate":  # 无机功能母体整名（P 中心无碳词干）：由 numbered 母体计数 + o_side 臂组装。
        from namepredict.layer5.phosphate import phosphate_names

        return phosphate_names(numbered)
    if kind in _EXO_SUF:  # 环外主基走 exocyclic worker（苯单取代→None 回落 chain_engine 保留名 variant）；开链同类无关（relation in_skeleton）。
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
    """词干以数字（1,3-thiazole）或 1H-（1H-pyrrole）开头时，与前缀需连字符分隔。"""
    return bool(body) and (body[0].isdigit() or body.startswith("1H-"))


def join_parent_name(prefix: str, parent: str) -> str:
    """拼接前缀与母体名（数字/1H- 前导时加连字符）。"""
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if _needs_join_hyphen(stem) else f"{prefix}{stem}"
    return f"{stereo}{body}"

def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str] | None:
    """拼接酯名：O 侧烷基作前缀、酸侧作主体（XX 酸 YY 酯）；无 O 侧取代基时输出 bare 酸酯名（benzoate/苯甲酸酯）。"""
    en, zh = names
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    if o:
        alk_en, alk_zh = o[0].get("en") or "", (o[0].get("zh") or "").rstrip("基")
    else:
        alk_en, alk_zh = "", ""
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


def _with_hydro_prefix(names: tuple[str, str], numbered: dict) -> tuple[str, str]:
    """把动态指示氢与 hydro 前缀依次拼到母体名前：顺序为 hydro + 指示氢 + 母体名（P-31.2.2，如 2,3-dihydro-1H-indole）；母体名已带静态 1H-（保留名）时不重复。"""
    parent = numbered.get("parent") or {}
    pre = parent.get("hydro_prefix")
    if not pre or not pre[0]:
        return names  # 仅氢化衍生物注入动态指示氢：全芳香环系的指示氢由保留名词干承载（1H-吡咯/9H-咔唑）或 fused_tree 路径（_ensure_fused_stem），此处注入会误产（[nH] 互变异构型得到 1H-pyridine）
    en, zh = names
    ind = parent.get("indicated_h") or ""
    if ind and not en.startswith(ind):
        en, zh = ind + en, ind + zh
    return (join_parent_name(pre[0], en), join_parent_name(pre[1], zh))


def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    """组装入口：取名 → 前缀 → 阴离子/R-S/金属盐后缀。"""
    from namepredict.layer5.stereo import apply_rs_prefix
    kind, n = _parent_n(numbered)
    if not _ensure_fused_stem(numbered):
        return _unsupported(n, kind)
    names = _names_for(kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    names = _with_hydro_prefix(names, numbered)

    joined = join_kind_name(kind, _prefix_for(numbered, kind, n), names, numbered)
    if joined is None:
        return _unsupported(n, kind)
    en, zh = joined
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
