"""L5 链式词干命名引擎：按 _Chain 规格数据驱动生成链状母体名（饱和/烯/炔/环）。"""
from __future__ import annotations
from dataclasses import KW_ONLY, dataclass, replace
from namepredict.layer5.stems import (
    ALKANE_EN, ALKANE_ZH, _en_stem, alkane_zh, zh_stem,
)
from namepredict.layer5.stereo import _ez_prefix, ez_for_parent
from namepredict.constants import MULT_EN, MULT_ZH, HALIDE_EN, HALO_ZH

def _pair(en_map: dict, zh_map: dict, n: int) -> tuple[str, str] | None:
    """从双语映射表取 n 的 (en, zh) 对，缺失返回 None。"""
    en, zh = en_map.get(n), zh_map.get(n); return (en, zh) if en and zh else None


def _alkane_names(n: int) -> tuple[str, str] | None:
    """查 C1–C35 烷烃英文/中文名表。"""
    return _pair(ALKANE_EN, ALKANE_ZH, n)


def _omit_term_locant(n: int, loc: int | None, omit: bool) -> bool:
    """端位/默认位次省略判定：omit、无位次或 C1–C2 的 1 位。"""
    return omit or loc is None or (loc == 1 and n <= 2)


def _NO_OMIT(n: int, loc: int | None, omit: bool) -> bool:
    """恒不省略位次的占位规则。"""
    return False


def _pair_loc_str(locs: list[int]) -> str:
    """位次列表拼接为逗号分隔字符串。"""
    return ",".join(str(x) for x in locs)


# C1/C2 开链英文 IUPAC 保留名（formic/acetic…）；C3+ 系统名由词干生成（_chain_plain 回落），中文无保留名。
_RETAINED = {
    "acid": {1: ("formic acid", "甲酸"), 2: ("acetic acid", "乙酸")},
    "acyl": {1: ("formyl", "甲酰基"), 2: ("acetyl", "乙酰基")},
    "aldehyde": {1: ("formaldehyde", "甲醛"), 2: ("acetaldehyde", "乙醛")},
    "amide": {1: ("formamide", "甲酰胺"), 2: ("acetamide", "乙酰胺")},
    "nitrile": {1: ("formonitrile", "甲腈"), 2: ("acetonitrile", "乙腈")},
    "ester": {1: ("formate", "甲酸"), 2: ("acetate", "乙酸")},
}


def _retained_plain(kind: str):
    """C1/C2 保留名 plain_fn：命中返回 (en,zh)，否则 None 回落词干生成。"""
    table = _RETAINED[kind]
    return lambda n: table.get(n)


def _ene_loc_kept(numbered: dict) -> int | None:
    """取应保留的烯位次：omit_ene_locant 开启时为 None。"""
    return None if numbered.get("omit_ene_locant") else numbered.get("ene_locant")


def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    """若母体含 E/Z 则给名称对加前缀。"""
    if pair is None: return None
    ez = ez_for_parent(numbered)
    return f"{ez}{pair[0]}", f"{ez}{pair[1]}"


def _fg_record(numbered: dict, kind: str) -> dict | None:
    """按 kind 查找主官能团位次记录（稀疏 fg_locants 列表）。"""
    return next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == kind), None)


def _parent_multiplicity(numbered: dict) -> int | None:
    """主官能团数量：facts.multiplicity，否则旧 principal_group_count。"""
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

def _has_ene(numbered: dict) -> bool:
    """判断母体是否存在烯键（含父字典 double_bond(s)）。"""
    p = numbered.get("parent") or {}
    return bool(
        numbered.get("ene_locant") or numbered.get("ene_locants")
        or p.get("double_bond") or p.get("double_bonds")
    )


def _has_yne(numbered: dict) -> bool:
    """判断母体是否存在炔键（yne_locant 或父字典 triple_bond）。"""
    p = numbered.get("parent") or {}
    return bool(numbered.get("yne_locant") or p.get("triple_bond"))


# ===== 链式词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环 (替代 if-kind 枚举) =====
def _chain_unsat(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """通用不饱和段引擎: 混合烯炔先组合段, 否则炔段优先、烯段其次 — 段式/融合式由 spec 数据驱动."""
    has_y = _has_yne(numbered) and (spec.yne_seg is not None or spec.yne_suf is not None)
    has_e = _has_ene(numbered) and (spec.ene_seg is not None or spec.ene_base is not None)
    if has_e and has_y:
        top = _chain_enyne(spec, n, numbered)
        if top is not None:
            return top
    if has_y:
        top = _chain_yne(spec, n, numbered)
        if top is not None:
            return top
    if has_e:
        top = _chain_ene(spec, n, numbered)
        if top is not None:
            return top
    return None

def _chain_enyne(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """混合烯炔段引擎：烯(en)前、炔(yne)后，合并两组位次；任一段为多键时词干加 "a"（P-31.1.1.2 euphonic a），烯段尾 e 在 yne 前省略。"""
    s, zs = _chain_stem_pair(spec, n)
    if s is None or zs is None:
        return None
    ene = numbered.get("ene_locant")
    enes = numbered.get("ene_locants")
    yne = numbered.get("yne_locant")
    ynes = numbered.get("yne_locants")
    ne = len(enes) if enes else (1 if ene is not None else 0)
    ny = len(ynes) if ynes else (1 if yne is not None else 0)
    if ne == 0 or ny == 0:
        return None
    multi = ne >= 2 or ny >= 2
    a_en = "a" if multi else ""
    e_loc = _pair_loc_str(enes) if ne >= 2 else str(ene)
    y_loc = _pair_loc_str(ynes) if ny >= 2 else str(yne)
    me, mz = MULT_EN.get(ne), MULT_ZH.get(ne)
    my_e, my_z = MULT_EN.get(ny), MULT_ZH.get(ny)
    if me is None or mz is None or my_e is None or my_z is None:
        return None
    if spec.ene_base is not None and spec.yne_suf is not None:   # 融合式: 烯段 MULT+en, 炔段 MULT+yne_suf
        e_seg = (f"{me}en", f"{mz}烯")
        y_suf = (f"{my_e}{spec.yne_suf[0]}", f"{my_z}{spec.yne_suf[1]}")
        return (f"{s}{a_en}-{e_loc}-{e_seg[0]}-{y_loc}-{y_suf[0]}",
                f"{zs}-{e_loc}-{e_seg[1]}-{y_loc}-{y_suf[1]}")
    if spec.unsat_polyol:  # 多 FG 词干模式: 混合烯段 elide e (pent-1-en-4-yne), FG 后缀由 _chain_names 拼接
        e_seg = (f"{me}en", f"{mz}烯")
        y_seg = (f"{my_e}yne", f"{my_z}炔")
        return (f"{s}{a_en}-{e_loc}-{e_seg[0]}-{y_loc}-{y_seg[0]}",
                f"{zs}-{e_loc}-{e_seg[1]}-{y_loc}-{y_seg[1]}")
    if spec.fg is not None:                                      # 段式: 烯段/炔段 + FG 位次
        fg = _fg_locant(numbered, spec.fg)
        if fg is None:
            return None
        e_seg = (f"{me}{spec.ene_seg[0]}", f"{mz}{spec.ene_seg[1]}")
        y_seg = (f"{my_e}{spec.yne_seg[0]}", f"{my_z}{spec.yne_seg[1]}")
        return (f"{s}{a_en}-{e_loc}-{e_seg[0]}-{y_loc}-{y_seg[0]}-{fg}-{spec.en_suf}",
                f"{zs}-{e_loc}-{e_seg[1]}-{y_loc}-{y_seg[1]}-{fg}-{spec.zh_suf}")
    return None

def _chain_yne(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """炔段引擎：融合式加炔后缀，段式（醇/酮）需 FG 位次；多炔用 MULT 后缀（hexa-1,5-diyne）。"""
    s, zs = _chain_stem_pair(spec, n)
    if s is None or zs is None:
        return None
    yne = numbered.get("yne_locant")
    ynes = numbered.get("yne_locants")
    multi = bool(ynes and len(ynes) >= 2)
    if spec.yne_suf is not None:      # 融合式: 炔后缀 (酸), omit 时无位次
        if n < 2:
            return None
        if multi:
            m_en, m_zh = MULT_EN.get(len(ynes)), MULT_ZH.get(len(ynes))
            if not m_en or not m_zh:
                return None
            loc = _pair_loc_str(ynes)
            suf = (f"{m_en}{spec.yne_suf[0]}", f"{m_zh}{spec.yne_suf[1]}")
            return f"{s}a-{loc}-{suf[0]}", f"{zs}-{loc}-{suf[1]}"
        omit = numbered.get("omit_yne_locant", False)
        if omit or yne is None:
            return f"{s}{spec.yne_suf[0]}", f"{zs}{spec.yne_suf[1]}"
        return f"{s}-{yne}-{spec.yne_suf[0]}", f"{zs}-{yne}-{spec.yne_suf[1]}"
    if spec.unsat_polyol:  # 多 FG 词干模式: 炔段保留 e (but-2-yne), FG 后缀由 _chain_names 拼接
        if yne is None and not ynes:
            return None
        if multi:
            m_en, m_zh = MULT_EN.get(len(ynes)), MULT_ZH.get(len(ynes))
            if not m_en or not m_zh:
                return None
            loc = _pair_loc_str(ynes)
            return (f"{s}a-{loc}-{m_en}yne", f"{zs}-{loc}-{m_zh}炔")
        if numbered.get("omit_yne_locant", False):
            return f"{s}yne", f"{zs}炔"
        return f"{s}-{yne}-yne", f"{zs}-{yne}-炔"
    fg = _fg_locant(numbered, spec.fg) if spec.fg else None   # 段式: 需 FG 位次 (醇/酮)
    if yne is None or fg is None:
        return None
    if multi:
        m_en, m_zh = MULT_EN.get(len(ynes)), MULT_ZH.get(len(ynes))
        if not m_en or not m_zh:
            return None
        loc = _pair_loc_str(ynes)
        seg = (f"{m_en}{spec.yne_seg[0]}", f"{m_zh}{spec.yne_seg[1]}")
        return (f"{s}a-{loc}-{seg[0]}-{fg}-{spec.en_suf}",
                f"{zs}-{loc}-{seg[1]}-{fg}-{spec.zh_suf}")
    return (
        f"{s}-{yne}-{spec.yne_seg[0]}-{fg}-{spec.en_suf}",
        f"{zs}-{yne}-{spec.yne_seg[1]}-{fg}-{spec.zh_suf}",
    )

def _fused_ene_suf(spec: "_Chain", m: int) -> tuple[str, str] | None:
    """融合式烯后缀: MULT[m] + 单烯基座 (acid→enoic/dienoic…); m 超界返回 None."""
    if m < spec.ene_m_min:
        return None
    me, mz = MULT_EN.get(m), MULT_ZH.get(m)
    if me is None or mz is None:
        return None
    return f"{me}{spec.ene_base[0]}", f"{mz}{spec.ene_base[1]}"

def _chain_ene(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """烯段引擎：多烯/单烯、融合式/段式按 spec 字段分支。"""
    s, zs = _chain_stem_pair(spec, n)
    if s is None or zs is None:
        return None
    enes = numbered.get("ene_locants")
    if enes and len(enes) >= 2:       # 多烯: 词干加 "a", 后缀生成式
        if spec.ene_base is not None:        # 融合式 (酸/醛/腈/二酸/酰胺/多烯)
            fused = _fused_ene_suf(spec, len(enes))
            if fused is None or n < spec.ene_n_min:
                return None
            ez = spec.ez_ene_multi(numbered) if spec.ez_ene_multi else ""
            loc = ",".join(str(x) for x in enes)
            return f"{ez}{s}a-{loc}-{fused[0]}", f"{ez}{zs}-{loc}-{fused[1]}"
        if spec.unsat_polyol:  # 多 FG 词干模式: 多烯 MULT+ene (buta-1,3-diene), FG 后缀由 _chain_names 拼接
            m_en, m_zh = MULT_EN.get(len(enes)), MULT_ZH.get(len(enes))
            if not m_en or not m_zh:
                return None
            ez = spec.ez_ene_multi(numbered) if spec.ez_ene_multi else ""
            loc = ",".join(str(x) for x in enes)
            return f"{ez}{s}a-{loc}-{m_en}ene", f"{ez}{zs}-{loc}-{m_zh}烯"
        if spec.fg is not None:              # 段式多烯 (醇/酮/胺/硫醇): MULT[m]+烯段
            fg = _fg_locant(numbered, spec.fg) if spec.fg else None
            if fg is None:
                return None
            m_en, m_zh = MULT_EN.get(len(enes)), MULT_ZH.get(len(enes))
            if not m_en or not m_zh:
                return None
            me, mz = f"{m_en}{spec.ene_seg[0]}", f"{m_zh}{spec.ene_seg[1]}"
            ez = spec.ez_ene_multi(numbered) if spec.ez_ene_multi else ""
            loc = ",".join(str(x) for x in enes)
            return f"{ez}{s}a-{loc}-{me}-{fg}-{spec.en_suf}", f"{ez}{zs}-{loc}-{mz}-{fg}-{spec.zh_suf}"
        return None
    ene = numbered.get("ene_locant")
    if spec.ene_base is not None:            # 融合式单烯 (酸/醛/腈/二酸/酰胺/多烯)
        if spec.ene_special is not None:
            sp = spec.ene_special(n, numbered)
            if sp is not None:
                return sp
        if ene is None or n < spec.ene_single_min:
            return None
        fused = _fused_ene_suf(spec, 1)
        if fused is None:
            return None
        ez = spec.ez_ene(numbered) if spec.ez_ene else ""
        if spec.ene_loc_omit and numbered.get("omit_ene_locant"):  # 环状单烯 (cyclohexene)：位次 1 隐含省略 (P-31.1.2)。
            return f"{ez}{s}{fused[0]}", f"{ez}{zs}{fused[1]}"
        return f"{ez}{s}-{ene}-{fused[0]}", f"{ez}{zs}-{ene}-{fused[1]}"
    if spec.ene_omit_aware:
        ene = _ene_loc_kept(numbered)
    if spec.unsat_polyol:  # 多 FG 词干模式: 烯段保留 e (but-2-ene), FG 后缀由 _chain_names 拼接
        if ene is None:
            return None
        ez = (spec.ez_ene(numbered) if spec.ez_ene else "") or ""
        return (f"{ez}{s}-{ene}-ene", f"{ez}{zs}-{ene}-烯")
    fg = _fg_locant(numbered, spec.fg) if spec.fg else None   # 段式单烯 (醇/酮)
    if ene is None or fg is None:
        return None
    ez = (spec.ez_ene(numbered) if spec.ez_ene else "") or ""
    if n <= 2 and ene == 1 and fg == 1:  # 短链 (C≤2) 单烯单 FG：烯只能 1(-2)、后缀锚定 1，位次无歧义省略并融合 (P-14.3.4.2/4.4；eth-1-en-1-amine → ethenamine，乙-1-烯-1-胺 → 乙烯胺)；另一端取代基（如 2-nitro）位次照常由前缀保留。
        return (
            f"{ez}{s}{spec.ene_seg[0]}{spec.en_suf}",
            f"{ez}{zs}{spec.ene_seg[1]}{spec.zh_suf}",
        )
    return (
        f"{ez}{s}-{ene}-{spec.ene_seg[0]}-{fg}-{spec.en_suf}",
        f"{ez}{zs}-{ene}-{spec.ene_seg[1]}-{fg}-{spec.zh_suf}",
    )

@dataclass(frozen=True)
class _Chain:
    """链式词干引擎配置规格：词缀后缀、位次规则、烯/炔段与 scaffold 变体等字段。"""
    kind: str
    en_suf: str            # 词缀后缀: "ol" / "one" / "ene" / "oic acid"
    zh_suf: str            # "醇" / "酮" / "烯" / "酸"
    coda: str = "an"       # 饱和词干后接 "an"; alkane/thiol/anhydride 特例 ""/"ane"
    _: KW_ONLY
    no_loc: str = "plain"  # "plain"=无位次输出普通名; "none"=返回 None (ketone)
    omit_rule: object = _NO_OMIT  # (n, loc, omit) -> bool  True=省略
    fg: str | None = None              # FG 类别名 (fg_locants 记录 kind): oh/amine/ketone/sh
    need: int | None = None            # FG 位次数要求 (need=1 单 FG; 多 FG 由数量生成)
    plain_maps: tuple | None = None    # 俗名表 (en_dict, zh_dict), 查不到时回落生成
    plain_fn: object = None            # 派生命名 (n)->pair — 酸酐从酸派生
    ene_seg: tuple = ("en", "烯")       # 段式单烯段; thiol 用 ("ene","烯") (P-57 保留 e)
    yne_seg: tuple = ("yn", "炔")       # 段式炔段; thiol 用 ("yne","炔")
    ene_base: tuple | None = None       # 融合式烯基座 = 单烯后缀 ("enoic acid","烯酸") — 实际后缀 = MULT[m]+基座
    ene_m_min: int = 1                  # 融合式烯最小烯数 (polyene/cyclopolyene 为 2)
    ene_special: object = None          # 融合式单烯俗名钩子 (n, numbered)->pair — amide 丙烯酰胺
    yne_suf: tuple | None = None        # 融合式炔后缀 ("ynoic acid","炔酸") — 酸
    ez_ene: object = None               # (numbered)->str  单烯 E/Z
    ez_ene_multi: object = None         # (numbered)->str  多烯 E/Z
    ene_n_min: int = 4                  # 多烯融合式 n 下限 (polyene 类为 0)
    ene_single_min: int = 2             # 融合式单烯 n 下限 (diacid 为 3)
    ene_omit_aware: bool = False        # 烯段受 omit_ene_locant 影响 (环系 FG)
    ene_loc_omit: bool = False          # 融合式单烯省略位次 (环单烯: cyclohexene)
    cyclic: bool = False                # 恒加环前缀 (纯烃环/环系 FG)
    cyclic_unsat: bool = False          # 仅烯/炔段时加环 (cyclopolyene: 无烯回落纯烷烃)
    zh_full: bool = False               # 中文词干保留完整烷烃后缀 "烷" (环烷/回落)
    wrap: object = None                 # (pair, numbered)->pair  整体包裹 (E/Z)
    unsat_polyol: bool = False          # 多 FG 词干支持烯/炔插入 (diol/triol: but-2-ene-1,4-diol)
    variant: dict | None = None         # {scaffold: {multiplicity: 生成式之上的特例字段覆盖}}; None 键=开链
    # assembler._names_for 按 scaffold_id 注入后一维化; 苯环保留名 (phenol/benzoic…)、acid 草酸特例。
    stem: tuple | None = None           # (en_stem, zh_stem) — 稠环/杂环 scaffold 词干覆盖 (naphthalen/萘…)
    mult_ok: bool = False               # 支持数量后缀生成 (acid/alcohol/amine/ketone)
    mult_zh_full: bool = False          # 多 FG 中文词干保留完整 "烷" (醇/胺)
    mult_unsat_polyol: bool = False     # 多 FG 词干支持烯/炔插入 (醇/胺/硫醇: but-2-ene-1,4-diol)
    aromatic: bool = False              # 芳香环 scaffold 标记 (由 assembler 注入); 醇→酚 语义在此消费

def _chain_zh_base(n: int) -> str | None:
    """中文烷烃全名去后缀得词干（丁烷→丁）。"""
    z = alkane_zh(n)
    return zh_stem(z) if z else None


def _chain_stem_pair(spec: "_Chain", n: int) -> tuple[str, str] | None:
    """取双语词干：unsat_polyol（多 FG）中文仅当 zh_full（醇/胺/硫醇）或环状时保留完整 "烷"，否则用去"烷"词干（戊-2,4-二酮）。"""
    if spec.stem:
        return spec.stem
    s = _en_stem(n)
    zs = alkane_zh(n) if (spec.unsat_polyol and (spec.zh_full or spec.cyclic)) else _chain_zh_base(n)
    if s is None or zs is None:
        return None
    return s, zs

def _elide_parent_e(stem: str, suffix: str) -> str:
    """母体氢化物结尾 'e' 的省略（P-60.2(a)）：仅后缀以元音 a/i/o/u/y 开头时省略，辅音结尾 (diol/dione/diamine/carbaldehyde) 保留（oxolane-3,4-diol）；判断下沉到拼接端、用最终 en_suf，避免 `_ring_stem` 无差别 rstrip('e') 的 bug。"""
    return stem[:-1] if stem.endswith("e") and suffix[:1].lower() in "aeiouy" else stem


def _chain_plain(spec: _Chain, s: str, zs: str, n: int) -> tuple[str, str]:
    """普通名: 俗名表 → 派生命名 → 词干拼接."""
    pair = _pair(*spec.plain_maps, n) if spec.plain_maps else None
    if pair is None and spec.plain_fn is not None:
        pair = spec.plain_fn(n)
    return pair if pair is not None else (f"{_elide_parent_e(s, spec.en_suf)}{spec.coda}{spec.en_suf}", f"{zs}{spec.zh_suf}")


def _radical_plain(n: int) -> tuple[str, str] | None:
    """烷基型省略形态：自由价在 1 位时 'ethan-1-yl'→'ethyl'（中文 '乙烷-1-基'→'乙基'，P-29.2 方法 1）；词干缺失回落 None 由调用方兜底。"""
    s, zs = _en_stem(n), _chain_zh_base(n)
    return (f"{s}yl", f"{zs}基") if s and zs else None

def _mult_elide(en_m: str, suf: str) -> str:
    """英文数量前缀元音省略：前缀尾 'a' 且后缀元音开头时去 'a'（P-14.3.2：tetra-ol→tetrol）；di/tri 不受影响。"""
    return en_m[:-1] if en_m.endswith("a") and suf[:1].lower() in "aeiou" else en_m


def _generated_mult_fields(spec: _Chain, mult: int) -> dict | None:
    """数量后缀生成: MULT[m] + 基础后缀; acid 特判烯基/炔基/俗名; 数量超 MULT 表返回 None."""
    en_m, zh_m = MULT_EN.get(mult), MULT_ZH.get(mult)
    if not en_m or not zh_m:
        return None
    fields = dict(
        en_suf=f"{_mult_elide(en_m, spec.en_suf)}{spec.en_suf}",
        zh_suf=f"{zh_m}{spec.zh_suf}",
        coda="ane", need=mult, no_loc="none",
        omit_rule=_NO_OMIT,
        plain_maps=None,
    )
    if spec.kind == "acid":  # 多酸烯基基座: enedioic/enetri…oic (保留 e); 中文 烯+数量酸。
        fields["ene_base"] = (f"ene{en_m}oic acid", f"烯{zh_m}酸")
        fields["yne_suf"] = None
        fields["ene_single_min"] = 3
    if spec.mult_zh_full:
        fields["zh_full"] = True
    if spec.mult_unsat_polyol:
        fields["unsat_polyol"] = True
    return fields


def _radical_terminal_yl_elide(pair: tuple[str, str], n: int) -> tuple[str, str]:
    """无环自由基自由价在 C-1 时省略 -1-（饱和链 P-29.2 方法 1 的烯/炔拓展）；丙-1-烯-1-基→丙-1-烯基、丁-3-烯-1-炔-1-基→丁-3-烯-1-炔基，C2 烯/炔（乙-1-炔-1-基→乙炔基、eth-1-yn-1-yl→ethynyl）位次 1 也省略，因短链位次无歧义（同 _omit_term_locant C1–C2 1 位）。"""
    en, zh = pair
    if en.endswith("-1-yl"):
        en = en[:-5] + "yl"
    if zh.endswith("-1-基"):
        zh = zh[:-4] + "基"
    if n <= 2:
        en = en.replace("-1-", "", 1)
        zh = zh.replace("-1-", "", 1)
    return en, zh


def _chain_names(spec: _Chain, n: int, numbered: dict) -> tuple[str, str] | None:
    """单链词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环; 烯/炔段插入由 spec 数据驱动."""
    if spec.aromatic and spec.kind == "alcohol":  # 芳香环醇统一"酚"（苯酚系；萘/吡啶/吲哚/喹啉同），主路径自动，非 variant 特例。
        spec = replace(spec, zh_suf="酚")
    mult = _parent_multiplicity(numbered)
    if mult is not None:
        if mult > 1 and spec.mult_ok:
            var = _generated_mult_fields(spec, mult)
            if var is None:
                return None
            var.update((spec.variant or {}).get(mult) or {})   # 特例覆盖 (acid 草酸/烯二酸)
            spec = replace(spec, **var)
        elif mult == 1:  # 单 FG 保留名覆盖 (苯环 → phenol/benzoic 等); 开链/无该 scaffold variant 时空。
            var = (spec.variant or {}).get(1)
            if var:
                spec = replace(spec, **var)
    alk = _alkane_names(n)
    if not alk:
        return None
    top = _chain_unsat(spec, n, numbered)
    if top is not None and spec.unsat_polyol:  # 多 FG 词干模式: 词干 + FG 位次 + 多 FG 后缀 (but-2-ene-1,4-diol)
        rec = _fg_record(numbered, spec.fg)
        locs = rec["locants"] if rec else None
        if not locs or (spec.need is not None and len(locs) != spec.need):
            top = None
        else:
            loc_s = ",".join(str(x) for x in locs)
            top = (f"{top[0]}-{loc_s}-{spec.en_suf}", f"{top[1]}-{loc_s}-{spec.zh_suf}")
    if top is not None and spec.kind == "radical" and spec.stem is None \
            and not spec.cyclic and not spec.cyclic_unsat:
        top = _radical_terminal_yl_elide(top, n)  # 开链自由基母体：自由价在 C-1 时省略 -1-（环自由基/稠环词干不受影响，见 _radical_terminal_yl_elide）。
    if top is not None:
        if spec.cyclic or spec.cyclic_unsat:
            top = (f"cyclo{top[0]}", f"环{top[1]}")
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
    if spec.fg is not None:      # FG 位次 (diol/dione/alcohol/...): 记录缺失或数量不符 → no_loc
        rec = _fg_record(numbered, spec.fg)
        locs = rec["locants"] if rec else None
        if not locs or (spec.need is not None and len(locs) != spec.need):
            if spec.no_loc == "none":
                return None
            pair = _chain_plain(spec, s, zs, n)
        elif spec.omit_rule(n, locs[0], rec["omit"] if rec else False):
            pair = _chain_plain(spec, s, zs, n)
        else:
            loc_s = ",".join(str(x) for x in locs)
            if spec.stem:  # 稠环 scaffold 词干已含完整基座（naphthalen/萘），直接拼后缀；结尾 'e' 按后缀首字母决定是否省略（P-60.2(a)，见 _elide_parent_e）。
                pair = (f"{_elide_parent_e(s, spec.en_suf)}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
            else:
                pair = (f"{s}{spec.coda}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
    else:
        pair = _chain_plain(spec, s, zs, n)
    if spec.cyclic:
        pair = (f"cyclo{pair[0]}", f"环{pair[1]}")
    return spec.wrap(pair, numbered) if spec.wrap is not None else pair


def _ac_hal_chain(hal_z: int) -> _Chain | None:
    """构造某卤素的酰卤链 spec（P-65.5）：后缀 -oyl halide/酰卤 与保留名（C1/C2、苯甲酰）都随卤素变。"""
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
                      None: {1: dict(plain_maps=None, plain_fn=lambda n, t=retained: t.get(n))},
                      "benzene": {1: dict(plain_maps=None,
                                          plain_fn=lambda n: (f"benzoyl {he}", f"苯甲酰{hz}"))},
                  })


# 按卤素原子序数索引的酰卤链 spec（F/Cl/Br/I）。
_ACYL_HALIDE_BY_HAL = {z: _ac_hal_chain(z) for z in HALIDE_EN}

_KIND_TABLE = {
    "alcohol": _Chain(kind="alcohol", en_suf="ol", zh_suf="醇",
                      fg="oh", need=1, omit_rule=_omit_term_locant,
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                      mult_ok=True, mult_zh_full=True, mult_unsat_polyol=True,
                      variant={
                          "benzene": {
                              1: dict(plain_maps=None, fg=None,
                                      plain_fn=lambda n: ("phenol", "苯酚")),
                          },
                      }),
    "ketone": _Chain(kind="ketone", en_suf="one", zh_suf="酮",
                     fg="ketone", need=1, no_loc="none",
                     omit_rule=lambda n, loc, omit: n <= 2 and loc == 1,
                     ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                     mult_ok=True,
                     mult_unsat_polyol=True),
    "alkane": _Chain(kind="alkane", en_suf="ane", zh_suf="烷", coda="",
                     omit_rule=lambda n, loc, omit: omit or n <= 3,
                     ene_base=("ene", "烯"), yne_suf=("yne", "炔"),
                     wrap=_with_ez),
    "acid": _Chain(kind="acid", en_suf="oic acid", zh_suf="酸",
                   ene_base=("enoic acid", "烯酸"),
                   yne_suf=("ynoic acid", "炔酸"),
                   ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent,
                   mult_ok=True,
                   variant={
                       None: {1: dict(plain_maps=None, plain_fn=_retained_plain("acid")),
                              2: dict(plain_maps=({2: "oxalic acid"}, {2: "草酸"}),
                                      yne_suf=None, ene_single_min=3)},
                       "benzene": {1: dict(plain_maps=None,
                                           plain_fn=lambda n: ("benzoic acid", "苯甲酸"))},
                   }),
    "ester": _Chain(kind="ester", en_suf="oate", zh_suf="酸",
                    ene_base=("enoate", "烯酸"),
                    yne_suf=("ynoate", "炔酸"),
                    ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent,
                    variant={
                        None: {1: dict(plain_maps=None, plain_fn=_retained_plain("ester"))},
                        "benzene": {1: dict(plain_maps=None,
                                            plain_fn=lambda n: ("benzoate", "苯甲酸"))},
                    }),
    # 酰基残基（P-65.1.7.2）：酸碳恒 locant 1、C3+ 系统名词干 coda "an" + "oyl"（propanoyl…），
    # C1/C2 走保留名（formyl/acetyl）；烯/炔与立体照 acid 融合式（enoyl/ynoyl）。
    "acyl": _Chain(kind="acyl", en_suf="oyl", zh_suf="酰基",
                   ene_base=("enoyl", "烯酰基"),
                   yne_suf=("ynoyl", "炔酰基"),
                   ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent,
                   variant={
                       None: {1: dict(plain_maps=None, plain_fn=_retained_plain("acyl"))},
                       # 苯环 exocyclic 酰基头 → 保留名 benzoyl/苯甲酰基（P-65.1.7.2 benzoic acid→benzoyl）；
                       # 杂环/碳环（furan-2-carbonyl/cyclopropanecarbonyl）由 assembler._exocyclic_acyl_names 出。
                       "benzene": {1: dict(plain_maps=None,
                                           plain_fn=lambda n: ("benzoyl", "苯甲酰基"))},
                   }),
    "thiol": _Chain(kind="thiol", en_suf="thiol", zh_suf="硫醇", coda="ane",
                    fg="sh", need=1, omit_rule=_omit_term_locant,
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                    mult_ok=True, mult_zh_full=True, mult_unsat_polyol=True,
                    ene_seg=("ene", "烯"), yne_seg=("yne", "炔")),
    "amine": _Chain(kind="amine", en_suf="amine", zh_suf="胺",
                    fg="amine", need=1, omit_rule=_omit_term_locant,
                    ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                    mult_ok=True, mult_zh_full=True, mult_unsat_polyol=True,
                    variant={
                        "benzene": {1: dict(plain_maps=None, fg=None,
                                            plain_fn=lambda n: ("aniline", "苯胺"))},
                    }),
    "aldehyde": _Chain(kind="aldehyde", en_suf="al", zh_suf="醛",
                       ene_base=("enal", "烯醛"),
                       yne_suf=("ynal", "炔醛"),
                       ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                       variant={
                           None: {1: dict(plain_maps=None, plain_fn=_retained_plain("aldehyde"))},
                           "benzene": {1: dict(plain_maps=None,
                                               plain_fn=lambda n: ("benzaldehyde", "苯甲醛"))},
                       }),
    "nitrile": _Chain(kind="nitrile", en_suf="enitrile", zh_suf="腈",
                      ene_base=("enenitrile", "烯腈"),
                      yne_suf=("ynenitrile", "炔腈"),
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent,
                      variant={
                          None: {1: dict(plain_maps=None, plain_fn=_retained_plain("nitrile"))},
                          "benzene": {1: dict(plain_maps=None,
                                              plain_fn=lambda n: ("benzonitrile", "苯甲腈"))},
                      }),
    "amide": _Chain(kind="amide", en_suf="amide", zh_suf="酰胺",
                    ene_base=("enamide", "烯酰胺"),

                    yne_suf=("ynamide", "炔酰胺"),
                    ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent,
                    variant={
                        None: {1: dict(plain_maps=None, plain_fn=_retained_plain("amide"))},
                        "benzene": {1: dict(plain_maps=None,
                                            plain_fn=lambda n: ("benzamide", "苯甲酰胺"))},
                    }),
    # 默认 chloride 由 _ac_hal_chain 生成，与按 parent.hal_z 选的 spec 同构（assembler._names_for 覆盖）。
    "acyl_halide": _ACYL_HALIDE_BY_HAL[17],
    "radical": _Chain(kind="radical", en_suf="yl", zh_suf="基", coda="an",
                      fg="radical", need=1, no_loc="none",
                      # 饱和无环链/单环烃自由价在 C-1 时省略位次（P-29.2 方法 1: ethyl/pentyl/2-phenylethyl）；
                      # 不饱和链（but-3-en-1-yl）与稠环/杂环（naphthalen-1-yl / pyridin-4-yl）走 unsat 段或 loc>1 保留。
                      omit_rule=lambda n, loc, omit: loc == 1,
                      plain_fn=_radical_plain,  # 省略位次用烷基型（ethyl/乙基），非 coda 拼接的 ethanyl/乙烷基
                      ene_seg=("en", "烯"), yne_seg=("yn", "炔"),
                      # 烯基自由基（取代基链含立体双键）需 E/Z 前缀：苯环母体上 prop-1-en-1-yl 等由
                      # 递归 * 锚定命名产出，切子分子已保立体（submol_build），此处按链位次拼 (1Z)-。
                      wrap=_with_ez,
                      variant={
                          "benzene": {1: dict(plain_maps=None,
                                              plain_fn=lambda n: ("phenyl", "苯基"),
                                              omit_rule=lambda n, loc, omit: True)},
                      }),
}
