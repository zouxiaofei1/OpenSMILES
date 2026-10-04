"""L5 链式词干命名引擎：按 _Chain 规格数据驱动生成链状母体名。"""
from __future__ import annotations
from dataclasses import KW_ONLY, dataclass, replace
from opensmiles.layer5.stems import _en_stem, alkane_en, alkane_zh, zh_stem
from opensmiles.layer5.stereo import _ez_prefix, _split_stereo_lead, ez_for_parent
from opensmiles.constants import CHAIN_RETAINED, Cl, EXO_RING_SUF, MULT_EN, MULT_ZH, HALIDE_EN, HALO_ZH, N, S

def _omit_term_locant(n: int, loc: int | None, omit: bool) -> bool:
    """端位/默认位次省略判定：omit、无位次或 C1–C2 的 1 位。"""
    return omit or loc is None or (loc == 1 and n <= 2)


def _NO_OMIT(n: int, loc: int | None, omit: bool) -> bool:
    """恒不省略位次的占位规则。"""
    return False


def _omit_all(n: int, loc: int | None, omit: bool) -> bool:
    """恒省略位次的占位规则。"""
    return True


def _const_plain(pair: tuple[str, str]):
    """构造忽略碳数的 plain_fn（保留名恒返回同一双语对）。"""
    return lambda n: pair


def _retained_plain(kind: str):
    """C1/C2 保留名 plain_fn：命中返回 (en,zh)。"""
    table = CHAIN_RETAINED[kind]
    return lambda n: table.get(n)


def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    """若母体含 E/Z 则给名称对加前缀。"""
    if pair is None: return None
    ez = ez_for_parent(numbered)
    return f"{ez}{pair[0]}", f"{ez}{pair[1]}"


def _fg_record(numbered: dict, kind: str) -> dict | None:
    """按 kind 查找主官能团位次记录（稀疏 fg_locants 列表）。"""
    return next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == kind), None)


def _fg_locs(numbered: dict, spec) -> tuple[dict | None, list | None]:
    """FG 位次记录与位次；缺失或不符时 locs 为 None。"""
    rec = _fg_record(numbered, spec.fg)
    locs = rec["locants"] if rec else None
    if not locs or (spec.need is not None and len(locs) != spec.need):
        return rec, None
    return rec, locs


def _parent_multiplicity(numbered: dict) -> int | None:
    """主官能团数量：facts.multiplicity 优先。"""
    parent = numbered.get("parent") or {}
    facts = parent.get("principal_expression_facts")
    if facts is not None:
        return facts.multiplicity
    count = parent.get("principal_group_count")
    return int(count) if count is not None else None


def _fg_locant(numbered: dict, kind: str) -> int | None:
    """单官能团位次（段式烯/炔、链名）；仅当恰好一个时返回，否则 None。"""
    rec = _fg_record(numbered, kind)
    locs = rec.get("locants") if rec else None
    return locs[0] if locs and len(locs) == 1 else None

def _yl_loc_omitted(spec: "_Chain", fg: int) -> bool:
    """自由价位次是否省略：无环链式母体的 C-1 自由价（P-29.2）。"""
    return (spec.yl_loc_omit and fg == 1
            and spec.stem is None and not spec.cyclic and not spec.cyclic_unsat)


def _fg_yl_tail(spec: "_Chain", fg: int) -> tuple[str, str]:
    """FG/自由价位次 + 后缀尾段（省略时位次不入串）。"""
    if _yl_loc_omitted(spec, fg):
        return spec.en_suf, spec.zh_suf
    return f"-{fg}-{spec.en_suf}", f"-{fg}-{spec.zh_suf}"


def _bond_locs(numbered: dict, b: str) -> tuple[int | None, list[int] | None, int, bool]:
    """取某类不饱和键的位次（首位次、列表、键数）。"""
    locs = numbered.get(f"{b}_locants") or None
    cnt = len(locs) if locs else 0
    return (locs[0] if locs else None), locs, cnt, bool(locs and len(locs) >= 2)


def _bond_form(spec: "_Chain", b: str) -> str | None:
    """不饱和段的形态：fused / polyol / seg 之一。"""
    if (spec.ene_base if b == "ene" else spec.yne_suf) is not None:
        return "fused"
    if spec.unsat_polyol:
        return "polyol"
    if spec.fg is not None:
        return "seg"
    return None


def _unsat_seg(spec: "_Chain", b: str) -> tuple:
    """段式烯/炔段：后缀首字母为辅音时保留末端 e（P-16.7.1(a)）。"""
    explicit = spec.ene_seg if b == "ene" else spec.yne_seg
    if explicit is not None:
        return explicit
    if spec.en_suf[:1].lower() not in "aeiouy":  # carboxylic acid / thiol 等保留 e
        return ("ene", "烯") if b == "ene" else ("yne", "炔")
    return ("en", "烯") if b == "ene" else ("yn", "炔")


def _bond_seg_str(spec: "_Chain", b: str, form: str, cnt: int, *, terminal: bool) -> tuple[str, str] | None:
    """不饱和键段的词段（不含词干与位次）。"""
    m_en, m_zh = MULT_EN.get(cnt), MULT_ZH.get(cnt)
    if m_en is None or m_zh is None:
        return None
    if form == "polyol":  # 词干内嵌: 烯/炔段保留 e，FG 后缀后拼
        if not terminal:
            return f"{m_en}en", f"{m_zh}烯"  # 连接段同融合式, 尾 e 在炔段前省略
        return (f"{m_en}ene", f"{m_zh}烯") if b == "ene" else (f"{m_en}yne", f"{m_zh}炔")
    if form == "fused":  # 融合式: 末段产出后缀基座
        if not terminal:
            return f"{m_en}en", f"{m_zh}烯"
        suf = spec.ene_base if b == "ene" else spec.yne_suf
        return (f"{m_en}{suf[0]}", f"{m_zh}{suf[1]}") if suf is not None else None
    seg = _unsat_seg(spec, b)
    return f"{m_en}{seg[0]}", f"{m_zh}{seg[1]}"


def _unsat_loc_omit(spec: "_Chain", b: str, form: str, single: bool, numbered: dict) -> bool:
    """该段位次是否省略（仅单段，P-14.3.4.2(d)）。"""
    if not single or form not in ("fused", "polyol"):
        return False
    if b == "ene":
        return (spec.ene_loc_omit or form == "polyol") and bool(numbered.get("omit_ene_locant"))
    return bool(numbered.get("omit_yne_locant")) and (spec.yne_loc_omit or form == "polyol")


def _chain_enyne(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """通用不饱和段引擎"""
    s, zs = _chain_stem_pair(spec, n)
    if s is None or zs is None:
        return None
    segs = [g for g in (("ene",) + _bond_locs(numbered, "ene"), ("yne",) + _bond_locs(numbered, "yne")) if g[3]]
    if not segs:
        return None
    if len(segs) == 2:  # 混合: 两段须共用同一形态，否则回落
        form = ("fused" if spec.ene_base is not None and spec.yne_suf is not None
                else "polyol" if spec.unsat_polyol
                else "seg" if spec.fg is not None else None)
    else:
        form = _bond_form(spec, segs[0][0])
    if form is None:
        return None
    single = len(segs) == 1

    ene_seg = next((g for g in segs if g[0] == "ene"), None)
    ez_hook = (spec.ez_ene_multi if ene_seg[4] else spec.ez_ene) if ene_seg is not None else None
    ez = ez_hook(numbered) if ez_hook else ""
    fg = _fg_locant(numbered, spec.fg) if form == "seg" else None
    if form == "seg":
        if fg is None:
            return None
        if single:  # 短链融合 / 环单烯（P-31.1.2）
            b, loc, _, cnt, _ = segs[0]
            seg = _bond_seg_str(spec, b, form, cnt, terminal=True)
            if seg is None:
                return None
            if n <= 2 and loc == 1 and fg == 1 and (b == "ene" or _yl_loc_omitted(spec, fg)):
                return f"{ez}{s}{seg[0]}{spec.en_suf}", f"{ez}{zs}{seg[1]}{spec.zh_suf}"
            if b == "ene" and cnt == 1 and spec.cyclic and spec.ene_loc_omit and loc == 1 and fg == 1:  # 仅单烯（多烯位次/倍增词不可省）
                loc_zh = "" if spec.zh_loc_omit else f"-{loc}-"
                return (f"{ez}{s}{seg[0]}-{fg}-{spec.en_suf}",
                        f"{ez}{zs}{loc_zh}{seg[1]}-{fg}-{spec.zh_suf}")
        tail_en, tail_zh = _fg_yl_tail(spec, fg)  # 段式: FG/自由价位次与后缀尾段（C-1 省略下沉在 _fg_yl_tail）
    else:
        tail_en = tail_zh = ""
    if spec.stem and s.endswith("ane"):  # 注入的完整母体名去 ane 后直接承接不饱和段（P-14.3.4.2）
        s = s[:-3]
    a_en = "a" if any(g[3] >= 2 for g in segs) else ""
    parts_en: list[str] = []
    parts_zh: list[str] = []
    for i, (b, loc, locs, cnt, multi) in enumerate(segs):
        seg = _bond_seg_str(spec, b, form, cnt, terminal=(i == len(segs) - 1))
        if seg is None:
            return None
        loc_s = None if _unsat_loc_omit(spec, b, form, single, numbered) else (
            ",".join(str(x) for x in locs) if multi else str(loc))
        parts_en.append(f"-{loc_s}-{seg[0]}" if loc_s is not None else seg[0])
        parts_zh.append(f"-{loc_s}-{seg[1]}" if loc_s is not None else seg[1])
    return (f"{ez}{s}{a_en}{''.join(parts_en)}{tail_en}",
            f"{ez}{zs}{''.join(parts_zh)}{tail_zh}")

@dataclass(frozen=True)
class _Chain:
    """链式词干引擎配置规格：词缀后缀、位次规则等字段。"""
    kind: str
    en_suf: str            # 词缀后缀: "ol" / "one" / "ene" / "oic acid"
    zh_suf: str            # "醇" / "酮" / "烯" / "酸"
    coda: str = "an"       # 饱和词干后接 "an"; 部分 kind 用 ""/"ane"
    _: KW_ONLY
    no_loc: str = "plain"  # "plain"=无位次普通名; "none"=返回 None
    omit_rule: object = _NO_OMIT  # (n, loc, omit) -> bool  True=省略
    fg: str | None = None              # FG 类别名 (fg_locants 记录 kind)
    need: int | None = None            # FG 位次数要求 (need=1 单 FG; 多 FG 由数量生成)
    plain_maps: tuple | None = None    # 俗名表 (en_dict, zh_dict), 查不到时回落生成
    plain_fn: object = None            # 派生命名 (n)->pair — 酸酐从酸派生
    ene_seg: tuple | None = None        # 段式单烯段; None = 按后缀首字母推导 (P-16.7.1(a))
    yne_seg: tuple | None = None        # 段式炔段; None 同上
    ene_base: tuple | None = None       # 融合式烯基座 = 单烯后缀 + 数量前缀
    yne_suf: tuple | None = None        # 融合式炔后缀 ("ynoic acid","炔酸") — 酸
    ez_ene: object = None               # (numbered)->str  单烯 E/Z
    ez_ene_multi: object = None         # (numbered)->str  多烯 E/Z
    ene_loc_omit: bool = False          # 融合式单烯省略位次 (乙烯 ethene / 环单烯 cyclohexene)
    zh_loc_omit: bool = True            # 中文省略不饱和位次 1；环外主基须保留
    yne_loc_omit: bool = False          # 融合式炔省略位次（P-14.3.4.2(d)）
    yl_loc_omit: bool = False           # 自由价在 C-1 时省略位次（P-29.2）
    cyclic: bool = False                # 恒加环前缀 (纯烃环/环系 FG)
    cyclic_unsat: bool = False          # 仅烯/炔段时加环 (cyclopolyene: 无烯回落纯烷烃)
    zh_full: bool = False               # 中文词干保留完整烷烃后缀 "烷" (多 FG: 丁烷-2,3-二醇; 回落)
    wrap: object = None                 # (pair, numbered)->pair  整体包裹 (E/Z)
    unsat_polyol: bool = False          # 多 FG 词干支持烯/炔插入 (diol/triol)
    variant: dict | None = None         # {scaffold: {multiplicity: 特例字段覆盖}}
    stem: tuple | None = None           # (en_stem, zh_stem) — 稠环/杂环词干覆盖
    aromatic: bool = False              # 芳香环 scaffold 标记；醇→酚 在此消费
    plain_hook: object = None           # (numbered)->pair|None  词尾由分子计数直接决定

def _chain_zh_base(n: int) -> str | None:
    """中文烷烃全名去后缀得词干（丁烷→丁）。"""
    z = alkane_zh(n)
    return zh_stem(z) if z else None


def _chain_stem_pair(spec: "_Chain", n: int) -> tuple[str, str] | None:
    """取双语词干：unsat_polyol 中文保留完整「烷」与否。"""
    if spec.stem:
        return spec.stem
    s = _en_stem(n)
    zs = alkane_zh(n) if (spec.unsat_polyol and (spec.zh_full or spec.cyclic)) else _chain_zh_base(n)
    if s is None or zs is None:
        return None
    return s, zs

def _elide_parent_e(stem: str, suffix: str) -> str:
    """母体结尾 'e' 省略（P-60.2(a)）：仅后缀元音开头时省。"""
    return stem[:-1] if stem.endswith("e") and suffix[:1].lower() in "aeiouy" else stem


def _chain_plain(spec: _Chain, s: str, zs: str, n: int) -> tuple[str, str]:
    """普通名: 俗名表 → 派生命名 → 词干拼接."""
    pair = None
    if spec.plain_maps: 
        en_map, zh_map = spec.plain_maps
        en, zh = en_map.get(n), zh_map.get(n)
        pair = (en, zh) if en and zh else None
    if pair is None and spec.plain_fn is not None:
        pair = spec.plain_fn(n)
    return pair if pair is not None else (f"{_elide_parent_e(s, spec.en_suf)}{spec.coda}{spec.en_suf}", f"{zs}{spec.zh_suf}")


def _hydrazide_plain(n: int) -> tuple[str, str] | None:
    """酰肼 C1/C2 保留名（P-66.3.1.2.1）。"""
    return {1: ("formohydrazide", "甲酰肼"), 2: ("acetohydrazide", "乙酰肼")}.get(n)


def hydrazide_chain_spec(spec: "_Chain") -> "_Chain":
    """酰胺 spec → 酰肼 spec（P-66.3.1.1）。"""
    return replace(spec, en_suf="hydrazide", zh_suf="酰肼", coda="ane",
                   ene_base=("enehydrazide", "烯酰肼"), yne_suf=("ynehydrazide", "炔酰肼"),
                   variant={None: {1: dict(plain_maps=None, plain_fn=_hydrazide_plain)}})


def _radical_plain(n: int) -> tuple[str, str] | None:
    """烷基型省略形态：自由价在 1 位时由词干拼 -yl/基（P-29.2）。"""
    s, zs = _en_stem(n), _chain_zh_base(n)
    return (f"{s}yl", f"{zs}基") if s and zs else None


def _generated_mult_fields(spec: _Chain, mult: int) -> dict | None:
    """数量后缀生成"""
    en_m, zh_m = MULT_EN.get(mult), MULT_ZH.get(mult)
    if not en_m or not zh_m:
        return None
    fields = dict(
        en_suf=f"{en_m[:-1] if en_m.endswith('a') and spec.en_suf[:1].lower() in 'aeiou' else en_m}{spec.en_suf}",  # P-14.3.2
        zh_suf=f"{zh_m}{spec.zh_suf}",
        coda="ane", need=mult, no_loc="none",
        omit_rule=_NO_OMIT,
        plain_maps=None,
    )
    if spec.kind == "acid":
        if spec.fg is not None and not spec.cyclic:  # 环外羧酸：后缀自带位次，不饱和段独立保留 e（P-65.1.1）
            fields.update(ene_base=None, yne_suf=None, unsat_polyol=True)
        else:  # 多酸烯基基座: 保留 e，中文 烯+数量酸
            fields["ene_base"] = (f"ene{en_m}oic acid", f"烯{zh_m}酸")
            fields["yne_suf"] = None
    if spec.kind == "ester":  # 多酯烯基/炔基基座: enedioate/ynedioate
        fields["ene_base"] = (f"ene{en_m}{spec.en_suf}", f"烯{zh_m}{spec.zh_suf}")
        fields["yne_suf"] = (f"yne{en_m}{spec.en_suf}", f"炔{zh_m}{spec.zh_suf}")
    return fields


def _free_valence_form(pair: tuple[str, str], order: int) -> tuple[str, str]:
    """自由价键级对应的基名形态（P-31.2.3）：双键亚基、三键次基。"""
    en_tail, zh_lead = ("ylidyne", "次") if order == 3 else ("ylidene", "亚")
    en, zh = pair
    if en.endswith("yl"):
        en = en[:-2] + en_tail
    if zh.endswith("基"):
        stem = zh[:-1]
        zh = f"{stem}{zh_lead}基" if stem.startswith("环") else f"{zh_lead}{stem}基"
    return en, zh


def _cyclo_stereo(pair: tuple[str, str]) -> tuple[str, str]:
    """环词干：前导 E/Z 立体块移到 cyclo/环 之前（P-91.2）。"""
    out: list[str] = []
    for text, cyc in zip(pair, ("cyclo", "环")):
        tag, rest = _split_stereo_lead(text)
        out.append(f"{tag}{cyc}{rest}")
    return out[0], out[1]


def _ring_prefix_located(numbered: dict) -> bool:
    """环上是否另带被编号前缀（P-66.6.1 主基位次不省）。"""
    chain = set((numbered.get("parent") or {}).get("chain") or [])
    return any(not s.get("o_side") and s.get("attach_idx") in chain
               for s in numbered.get("substituents") or [])


def _exo_ring_spec(spec: "_Chain", n: int, numbered: dict) -> "_Chain":
    """环外主基的 spec 改写（P-65.2.2 / P-66.6.1.1.3）。"""
    suf = EXO_RING_SUF.get(spec.kind)
    parent = numbered.get("parent") or {}
    if suf is not None and spec.kind == "ester" and parent.get("thio_side"):  # P-65.6.3.3.7.1 环外硫代羧酸 S-酯
        suf = (("carbothioate", "硫代甲酸"), None)
    if suf is not None and spec.kind == "amide" and parent.get("amide_z") == S:  # P-43 类 16：环外硫代酰胺
        suf = (("carbothioamide", "硫代甲酰胺"), ("carbothioamide", "硫代甲酰胺"))
    elif suf is not None and spec.kind == "amide" and parent.get("amide_z") == N:  # P-43 类 17：环外脒
        suf = (("carboximidamide", "亚氨酰胺"), ("carboximidamide", "亚氨酰胺"))
    elif suf is not None and spec.kind == "amide" and parent.get("hydrazide_n_idx") is not None:
        suf = (("carbohydrazide", "甲酰肼"), ("carbohydrazide", "甲酰肼"))  # P-66.3.1.1：环上 -CO-NHNH2 用 carbohydrazide
    facts = parent.get("principal_expression_facts")
    if suf is None or facts is None or facts.relation.value != "exocyclic":
        return spec
    singular, plural = suf
    if spec.kind == "acyl_halide":  # 卤素词随卤原子变化，环外后缀须动态拼接
        he, hz = HALIDE_EN.get(parent.get("hal_z")), HALO_ZH.get(parent.get("hal_z"))
        if he is None or hz is None:
            return spec
        singular = (f"{singular[0]} {he}", f"{singular[1]}{hz}")
    mult = facts.multiplicity
    if mult > 1 and plural is None:  # 酯/酰胺/腈/酰基头无多取代系统名。
        return spec
    sid = parent.get("scaffold_id")
    if mult == 1 and sid == "benzene":  # 苯单取代保留名由 variant 承担。
        return spec
    fields = dict(en_suf=singular[0], zh_suf=singular[1], fg=spec.kind, ene_base=None, yne_suf=None,
                 )  # 多取代后缀由数量机制生成
    if sid != "carbocycle" or parent.get("fused_tree") or parent.get("stem_generated"):  # 稠环/杂环/生成式词干已注入，位次恒显式。
        return replace(spec, **fields)
    base = (alkane_en(n), alkane_zh(n))
    if base[0] is None or base[1] is None:
        return spec
    if numbered.get("ene_locants") or numbered.get("yne_locants"):  # 环内不饱和：段式后缀，单烯 1 位 EN 省略
        fields.update(stem=None, coda="ane", cyclic=True, zh_loc_omit=False)
    else:  # 饱和环母体取完整氢化物（词干已含环前缀，故关 cyclic）。
        fields.update(stem=(f"cyclo{base[0]}", f"环{base[1]}"), coda="", cyclic=False)
    fields["omit_rule"] = lambda n, loc, omit: omit or not _ring_prefix_located(numbered)  # 干净环省略主基位次（cyclohexanecarboxylic acid）
    return replace(spec, **fields)


_ZH_CODA_KINDS = frozenset({"sulfonic", "sulfonate", "sulfonamide", "sulfonyl_chloride"})  # 中文带「烷」coda 的 S 链后缀母体


def _chain_names(spec: _Chain, n: int, numbered: dict) -> tuple[str, str] | None:
    """单链词干引擎: 数词干 + coda + 后缀 + 位次 + 环。"""

    if spec.plain_hook is not None:  # 无碳链母体（磷酸等 P 中心）：词尾由计数直接产出，不查碳数词表。
        hook = spec.plain_hook(numbered)
        if hook is not None:
            return hook
    if spec.aromatic and spec.kind == "alcohol":  # 芳香环醇统一"酚"，主路径自动。
        spec = replace(spec, zh_suf="酚")
    spec = _exo_ring_spec(spec, n, numbered)
    mult = _parent_multiplicity(numbered)
    if mult is not None and mult > 1:
        var = _generated_mult_fields(spec, mult)
        if var is None:
            return None
        var.update((spec.variant or {}).get(mult) or {})   # 特例覆盖 (acid 草酸/烯二酸)
        spec = replace(spec, **var)
    else:  # 非多 FG: 关闭多 FG 字段，再套用单 FG 保留名覆盖
        spec = replace(spec, unsat_polyol=False, zh_full=False)
        var = (spec.variant or {}).get(1) if mult == 1 else None
        if var:
            spec = replace(spec, **var)

    free_order = (numbered.get("parent") or {}).get("free_valence_order") or 0
    if spec.kind != "radical":  # 自由价形态只对自由基类母体生效
        free_order = 0
    top = _chain_enyne(spec, n, numbered)
    if top is not None and spec.unsat_polyol:  # 多 FG 词干模式: 词干 + FG 位次 + 后缀
        _, locs = _fg_locs(numbered, spec)
        if locs is None:
            top = None
        else:
            loc_s = ",".join(str(x) for x in locs)
            top = (f"{top[0]}-{loc_s}-{spec.en_suf}", f"{top[1]}-{loc_s}-{spec.zh_suf}")
    if top is not None:  # 无环自由基位次省略已下沉到段式引擎，环自由基不受影响。
        if spec.cyclic or spec.cyclic_unsat:
            top = _cyclo_stereo(top)
        if free_order > 1:
            top = _free_valence_form(top, free_order)
        return spec.wrap(top, numbered) if spec.wrap is not None else top
    if spec.stem:
        s, zs = spec.stem
    else:
        s = _en_stem(n)
        if spec.zh_full:
            zf = alkane_zh(n)
            zs = zf if zf else _chain_zh_base(n)
        else:
            zs = _chain_zh_base(n)
    if s is None or zs is None:
        return None
    if spec.fg is not None:      # FG 位次: 记录缺失或数量不符 → no_loc
        rec, locs = _fg_locs(numbered, spec)
        if locs is None:
            if spec.no_loc == "none":
                return None
            pair = _chain_plain(spec, s, zs, n)
        elif spec.omit_rule(n, locs[0], rec["omit"] if rec else False):
            pair = _chain_plain(spec, s, zs, n)
        else:
            loc_s = ",".join(str(x) for x in locs)
            if spec.stem:  # 稠环词干已含基座，直接拼后缀（P-60.2(a)）
                pair = (f"{_elide_parent_e(s, spec.en_suf)}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
            else:
                zstem = alkane_zh(n) if spec.kind in _ZH_CODA_KINDS else zs  # S 链后缀中文保留「烷」（丙烷-1-磺酰氯）
                pair = (f"{s}{spec.coda}-{loc_s}-{spec.en_suf}", f"{zstem}-{loc_s}-{spec.zh_suf}")
    else:
        pair = _chain_plain(spec, s, zs, n)
    if spec.cyclic:
        pair = (f"cyclo{pair[0]}", f"环{pair[1]}")
    if free_order > 1:
        pair = _free_valence_form(pair, free_order)
    return spec.wrap(pair, numbered) if spec.wrap is not None else pair


def _ac_hal_chain(hal_z: int) -> _Chain | None:
    """构造某卤素的酰卤链 spec（P-65.5）。"""
    he = HALIDE_EN.get(hal_z)
    hz = HALO_ZH.get(hal_z)
    if he is None or hz is None:
        return None
    retained = {1: (f"formyl {he}", f"甲酰{hz}"), 2: (f"acetyl {he}", f"乙酰{hz}")}
    return _Chain(kind="acyl_halide", en_suf=f"oyl {he}", zh_suf=f"酰{hz}",
                  ene_base=(f"enoyl {he}", f"烯酰{hz}"),
                  yne_suf=(f"ynoyl {he}", f"炔酰{hz}"),
                  ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent,
                  variant={
                      None: {1: dict(plain_maps=None, plain_fn=lambda n, t=retained: t.get(n)),
                             2: dict(plain_maps=None, plain_fn=_const_plain((f"oxalyl di{he}", f"草酰二{hz}")))},  # P-65.1.1 保留名：乙二酰二卤 = oxalyl
                      "benzene": {1: dict(plain_maps=None, plain_fn=_const_plain((f"benzoyl {he}", f"苯甲酰{hz}")))},
                  })

_ACYL_HALIDE_BY_HAL = {z: _ac_hal_chain(z) for z in HALIDE_EN}  # 按卤素原子序数索引的酰卤链 spec（F/Cl/Br/I）。

_OXO_TAIL = {  # P-66.1.1/P-67.1.3：oxo_kind × 氢数 → 功能母体词尾
    ("phosphate", 3): ("phosphoric acid", "磷酸"),
    ("phosphate", 2): ("dihydrogen phosphate", "磷酸二氢"),
    ("phosphate", 1): ("hydrogen phosphate", "磷酸氢"),
    ("phosphate", 0): ("phosphate", "磷酸"),
    ("phosphonate", 2): ("phosphonic acid", "膦酸"),
    ("phosphonate", 1): ("hydrogen phosphonate", "膦酸氢"),
    ("phosphonate", 0): ("phosphonate", "膦酸"),
    ("sulfate", 2): ("sulfuric acid", "硫酸"),
    ("sulfate", 1): ("hydrogen sulfate", "硫酸氢"),
    ("sulfate", 0): ("sulfate", "硫酸"),
    ("boronic", 2): ("boronic acid", "硼酸"),      # P-68.2.1：R-B(OH)2 以功能母体词尾命名，R 退为取代基
    ("boronic", 1): ("hydrogen boronate", "硼酸氢"),
    ("boronic", 0): ("boronate", "硼酸"),
}


def _oxoacid_tail(numbered: dict) -> tuple[str, str] | None:
    """含氧酸中心母体词尾：按 oxo_kind 与酸式氢数查表，表外返回 None。"""
    parent = numbered.get("parent") or {}
    return _OXO_TAIL.get((parent.get("oxo_kind"), int(parent.get("n_oh") or 0)))


def _benzene_retained(en: str, zh: str, *, fg_drop: bool = False, omit_all: bool = False) -> dict:
    """苯单取代保留名 variant（P-61.2）：C1 直接取保留名。"""
    fields = dict(plain_maps=None, plain_fn=_const_plain((en, zh)))
    if fg_drop:
        fields["fg"] = None
    if omit_all:
        fields["omit_rule"] = _omit_all
    return {1: fields}


_BENZENE_RETAINED = {  # 苯环单取代 scaffold 专属保留名（消费点 assembler）。
    "alcohol": _benzene_retained("phenol", "苯酚", fg_drop=True),
    "amine": _benzene_retained("aniline", "苯胺", fg_drop=True),
    "radical": _benzene_retained("phenyl", "苯基", omit_all=True),
    "acid": _benzene_retained("benzoic acid", "苯甲酸"),
    "sulfonic": _benzene_retained("benzenesulfonic acid", "苯磺酸", omit_all=True),
    "sulfonate": _benzene_retained("benzenesulfonate", "苯磺酸", omit_all=True),
    "sulfonamide": _benzene_retained("benzenesulfonamide", "苯磺酰胺", omit_all=True),
    "sulfonyl_chloride": _benzene_retained("benzenesulfonyl chloride", "苯磺酰氯", omit_all=True),
    "ester": _benzene_retained("benzoate", "苯甲酸"),
    "acyl": _benzene_retained("benzoyl", "苯甲酰基"),
    "aldehyde": _benzene_retained("benzaldehyde", "苯甲醛"),
    "nitrile": _benzene_retained("benzonitrile", "苯甲腈"),
    "amide": _benzene_retained("benzamide", "苯甲酰胺"),
}

_KIND_TABLE = {
    "alcohol": _Chain(kind="alcohol", en_suf="ol", zh_suf="醇",
                      fg="alcohol", need=1, omit_rule=_omit_term_locant,
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                      zh_full=True, unsat_polyol=True),
    "ketone": _Chain(kind="ketone", en_suf="one", zh_suf="酮",
                     fg="ketone", need=1, no_loc="none",
                     omit_rule=lambda n, loc, omit: n <= 2 and loc == 1,
                     ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,unsat_polyol=True),
    "thione": _Chain(kind="thione", en_suf="thione", zh_suf="硫酮", coda="ane",  # P-64.6.1：酮的硫族类似物，母体保留 e
                     fg="thione", need=1, omit_rule=_NO_OMIT,
                     ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent, unsat_polyol=True),
    "alkane": _Chain(kind="alkane", en_suf="ane", zh_suf="烷", coda="",
                     omit_rule=lambda n, loc, omit: omit or n <= 3,
                     ene_base=("ene", "烯"), yne_suf=("yne", "炔"),
                     ene_loc_omit=True, yne_loc_omit=True,
                     wrap=_with_ez),
    "acid": _Chain(kind="acid", en_suf="oic acid", zh_suf="酸",
                   ene_base=("enoic acid", "烯酸"),
                   yne_suf=("ynoic acid", "炔酸"),
                   ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                   variant={
                       None: {1: dict(plain_maps=None, plain_fn=_retained_plain("acid")),
                              2: dict(plain_maps=({2: "oxalic acid"}, {2: "草酸"}), yne_suf=None)},
                   }),
    "sulfonic": _Chain(kind="sulfonic", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                       en_suf="sulfonic acid", zh_suf="磺酸", coda="ane",
                       fg="oxoacid", need=1, omit_rule=_omit_term_locant),  # P-65.3.1 取代式：链/环母体 + 磺酸后缀
    "sulfonate": _Chain(kind="sulfonate", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                        en_suf="sulfonate", zh_suf="磺酸", coda="ane",
                        fg="oxoacid", need=1, omit_rule=_omit_term_locant),  # 磺酸酯：O 侧臂 + 磺酸酯后缀
    "sulfonamide": _Chain(kind="sulfonamide", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                          en_suf="sulfonamide", zh_suf="磺酰胺", coda="ane",
                          fg="sulfonamide", need=1, omit_rule=_omit_term_locant),
    "sulfonyl_chloride": _Chain(kind="sulfonyl_chloride", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                                en_suf="sulfonyl chloride", zh_suf="磺酰氯", coda="ane",
                                fg="oxoacid", need=1, omit_rule=_omit_term_locant),
    "ester": _Chain(kind="ester", en_suf="oate", zh_suf="酸",
                    ene_base=("enoate", "烯酸"),
                    yne_suf=("ynoate", "炔酸"),
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                    variant={
                        None: {1: dict(plain_maps=None, plain_fn=_retained_plain("ester")),
                               2: dict(plain_maps=({2: "oxalate"}, {2: "草酸"}))},  # P-65.1.1 保留名：乙二酸二酯 = oxalate
                    }),
    "phosphate": _Chain(kind="phosphate", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent, en_suf="phosphate", fg="oxoacid", zh_suf="磷酸", coda="",  # P-67.1.3 无机功能母体：词尾由 plain_hook 切换
                        plain_hook=_oxoacid_tail),
    "phosphonate": _Chain(kind="phosphonate", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent, en_suf="phosphonic acid", fg="oxoacid", zh_suf="膦酸", coda="",  # P-67.1.2：膦酸功能母体，碳臂作前缀
                          plain_hook=_oxoacid_tail),
    "sulfate": _Chain(kind="sulfate", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent, en_suf="sulfate", fg="oxoacid", zh_suf="硫酸", coda="",  # 硫酸酯功能母体
                      plain_hook=_oxoacid_tail),
    "boronic": _Chain(kind="boronic", ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent, en_suf="boronic acid", fg="oxoacid", zh_suf="硼酸", coda="",  # P-68.2.1：硼酸功能母体，碳基作前缀
                      plain_hook=_oxoacid_tail),
    "acyl": _Chain(kind="acyl", en_suf="oyl", zh_suf="酰基",  # 酰基残基（P-65.1.7.2）：酸碳恒 locant 1，C1/C2 走保留名
                   ene_base=("enoyl", "烯酰基"),
                   yne_suf=("ynoyl", "炔酰基"),
                   ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                   variant={None: {1: dict(plain_maps=None, plain_fn=_retained_plain("acyl"))}}),
    "thiol": _Chain(kind="thiol", en_suf="thiol", zh_suf="硫醇", coda="ane",
                    fg="thiol", need=1, omit_rule=_omit_term_locant,
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                     zh_full=True, unsat_polyol=True),
    "amine": _Chain(kind="amine", en_suf="amine", zh_suf="胺",
                    fg="amine", need=1, omit_rule=_omit_term_locant,
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                    zh_full=True, unsat_polyol=True),
    "aldehyde": _Chain(kind="aldehyde", en_suf="al",zh_suf="醛",
                       ene_base=("enal", "烯醛"),
                       yne_suf=("ynal", "炔醛"),
                       ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                       variant={None: {1: dict(plain_maps=None, plain_fn=_retained_plain("aldehyde"))}}),
    "nitrile": _Chain(kind="nitrile", en_suf="enitrile", zh_suf="腈",
                      ene_base=("enenitrile", "烯腈"),
                      yne_suf=("ynenitrile", "炔腈"),
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                      variant={None: {1: dict(plain_maps=None, plain_fn=_retained_plain("nitrile"))}}),
    "amide": _Chain(kind="amide", en_suf="amide", zh_suf="酰胺", 
                    ene_base=("enamide", "烯酰胺"),
                    yne_suf=("ynamide", "炔酰胺"),
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                    variant={
                        None: {1: dict(plain_maps=None, plain_fn=_retained_plain("amide")),
                               2: dict(plain_maps=({2: "oxamide"}, {2: "草酰胺"}))},  # P-66.1.1 保留名：乙二酰胺 = oxamide
                    }),
    "acyl_halide": _ACYL_HALIDE_BY_HAL[Cl],  # 默认 chloride spec，由 assembler 按 hal_z 覆盖
    "radical": _Chain(kind="radical", en_suf="yl", zh_suf="基", coda="an",
                      fg="radical", need=1, no_loc="none", yl_loc_omit=True,
                      omit_rule=lambda n, loc, omit: loc == 1,  # 饱和无环链/单环烃自由价在 C-1 时省略位次（P-29.2）
                      plain_fn=_radical_plain, 
                      wrap=_with_ez),  # 烯基自由基需 E/Z 前缀
}
