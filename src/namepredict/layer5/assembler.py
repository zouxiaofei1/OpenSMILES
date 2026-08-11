from __future__ import annotations
from dataclasses import dataclass
from namepredict.layer5.stems import (
    ACID_EN, ACID_ZH, ALCOHOL_EN, ALCOHOL_ZH,
    ALDEHYDE_EN, ALDEHYDE_ZH, ALKANE_EN, ALKANE_ZH,
    AMIDE_EN, AMIDE_ZH, ESTER_ACYL_EN,
    NITRILE_EN, NITRILE_ZH, SULFIDE_ALKYL_EN,
    SULFIDE_ALKYL_ZH, SULFIDE_SYM_EN, SULFIDE_SYM_ZH, _en_stem,
    alkane_zh, maybe_anion_names, maybe_metal_salt_names, zh_stem,
)
from namepredict.layer5.benzene_names import (
    benzene_parent_names, phenyl_parent_names,
    benzenediol_names,
    join_kind_name, pyridine_kind_names,
)
from namepredict.layer5.stereo_ez import _ez_prefix, ez_for_parent
from namepredict.layer5.unsat_acid import alkenamide_names
from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.types import NameResult
def _fail(meta: dict | None = None) -> NameResult: return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})
def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult: return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)
def _pair(en_map: dict, zh_map: dict, n: int) -> tuple[str, str] | None:
    en, zh = en_map.get(n), zh_map.get(n); return (en, zh) if en and zh else None
def _alkane_names(n: int) -> tuple[str, str] | None: return _pair(ALKANE_EN, ALKANE_ZH, n)
def _omit_term_locant(n: int, loc: int | None, omit: bool) -> bool:
    return omit or loc is None or (loc == 1 and n <= 2)
def _pair_loc_str(locs: list[int]) -> str:
    return ",".join(str(x) for x in locs)
def _acid_names(n: int) -> tuple[str, str] | None:
    return _pair(ACID_EN, ACID_ZH, n)
def _anhydride_from_acid(n: int) -> tuple[str, str] | None:
    plain = _acid_names(n)
    if not plain:
        return None
    en, zh = plain
    return en.replace(" acid", " anhydride"), f"{zh}酐"
def _ene_loc_kept(numbered: dict) -> int | None:
    return None if numbered.get("omit_ene_locant") else numbered.get("ene_locant")

def _sym_sulfide_names(n: int) -> tuple[str, str] | None:
    return _pair(SULFIDE_SYM_EN, SULFIDE_SYM_ZH, n)
def _asym_sulfide_names(n1: int, n2: int) -> tuple[str, str] | None:
    en1, en2 = SULFIDE_ALKYL_EN.get(n1), SULFIDE_ALKYL_EN.get(n2)
    zh1, zh2 = SULFIDE_ALKYL_ZH.get(n1), SULFIDE_ALKYL_ZH.get(n2)
    if not en1 or not en2 or not zh1 or not zh2:
        return None
    a, b = sorted([(en1, zh1), (en2, zh2)], key=lambda x: x[0])
    return f"{a[0]} {b[0]} sulfide", f"{a[1]}{b[1]}硫醚"
def _sulfide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    ns = parent.get("alkyl_ns")
    if not ns or len(ns) != 2:
        return None
    n1, n2 = int(ns[0]), int(ns[1])
    if n1 == n2:
        return _sym_sulfide_names(n1)
    return _asym_sulfide_names(n1, n2)

def _typed_acid_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "acid":
        return kind
    if facts.relation.value == "exocyclic" and kind == "cycloalkane":
        return "cycloalkanecarboxylic" if facts.multiplicity == 1 else "cycloalkane_polycarboxylic"
    if facts.relation.value == "exocyclic":
        return kind
    return "acid" if facts.multiplicity == 1 else "diacid" if facts.multiplicity == 2 else "polycarboxylic"

def _typed_ketone_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ketone":
        return kind
    if kind in {"cycloalkane", "cycloketone", "cycloalkanedione"}:
        return "cycloketone" if facts.multiplicity == 1 else "cycloalkanedione"
    if kind not in {"ketone", "dione"}:
        return kind
    return "ketone" if facts.multiplicity == 1 else "dione"

def _typed_ring_alcohol_kind(kind: str, numbered: dict, facts) -> str | None:
    if kind in {"cycloalkane", "cycloalcohol", "cycloalkanediol"}:
        return "cycloalcohol" if facts.multiplicity == 1 else "cycloalkanediol"
    if kind in {"benzene", "phenol", "benzenediol"}:
        return "phenol" if facts.multiplicity == 1 else "benzenediol"
    parent = numbered.get("parent") or {}
    scaffold = parent.get("scaffold_identity")
    if parent.get("typed_ring_expression_supported") and scaffold and scaffold.id == "naphthalene":
        return "naphthalenol" if facts.multiplicity == 1 else "naphthalenediol"
    return None

def _typed_alcohol_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "alcohol":
        return kind
    ring_kind = _typed_ring_alcohol_kind(kind, numbered, facts)
    if ring_kind:
        return ring_kind
    if kind not in {"alcohol", "diol", "triol"}:
        return kind
    return "alcohol" if facts.multiplicity == 1 else "diol" if facts.multiplicity == 2 else "triol"

def _typed_amine_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amine":
        return kind
    if kind in {"cycloalkane", "cycloamine"} and facts.multiplicity == 1:
        return "cycloamine"
    if kind not in {"amine", "diamine", "triamine", "tetraamine"}:
        return kind
    return {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"}.get(facts.multiplicity, kind)

def _typed_expression_kind(kind: str, numbered: dict) -> str:
    kind = _typed_acid_kind(kind, numbered)
    kind = _typed_ketone_kind(kind, numbered)
    kind = _typed_alcohol_kind(kind, numbered)
    return _typed_amine_kind(kind, numbered)
def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    if pair is None: return None
    from namepredict.layer5.stereo_ez import ez_for_parent
    ez = ez_for_parent(numbered)
    return f"{ez}{pair[0]}", f"{ez}{pair[1]}"

# ===== 链式词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环 (替代 if-kind 枚举) =====
def _chain_unsat(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """通用不饱和段引擎: 炔段优先, 烯段其次 — 段式(醇/酮)与融合式(酸)均由 spec 数据驱动."""
    from namepredict.layer5.unsat_acid import _has_ene, _has_yne
    if _has_yne(numbered) and (spec.yne_seg is not None or spec.yne_suf is not None):
        top = _chain_yne(spec, n, numbered)
        if top is not None:
            return top
    if _has_ene(numbered) and (spec.ene_seg is not None or spec.ene_base is not None):
        top = _chain_ene(spec, n, numbered)
        if top is not None:
            return top
    return None

def _chain_yne(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    s, zs = _en_stem(n), _chain_zh_base(n)
    if s is None or zs is None:
        return None
    yne = numbered.get("yne_locant")
    if spec.yne_suf is not None:      # 融合式: 炔后缀 (酸), omit 时无位次
        if n < 2:
            return None
        omit = numbered.get("omit_yne_locant", False)
        if omit or yne is None:
            return f"{s}{spec.yne_suf[0]}", f"{zs}{spec.yne_suf[1]}"
        return f"{s}-{yne}-{spec.yne_suf[0]}", f"{zs}-{yne}-{spec.yne_suf[1]}"
    fg = numbered.get(spec.loc) if spec.loc else None   # 段式: 需 FG 位次 (醇/酮)
    if yne is None or fg is None:
        return None
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
    s, zs = _en_stem(n), _chain_zh_base(n)
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
        if spec.mult_seg is not None:        # 段式多烯 (醇)
            fg = numbered.get(spec.loc) if spec.loc else None
            if fg is None:
                return None
            me, mz = spec.mult_seg.get(len(enes), ("", ""))
            if not me or not mz:
                return None
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
        return f"{ez}{s}-{ene}-{fused[0]}", f"{ez}{zs}-{ene}-{fused[1]}"
    if spec.ene_omit_aware:
        ene = _ene_loc_kept(numbered)
    fg = numbered.get(spec.loc) if spec.loc else None   # 段式单烯 (醇/酮)
    if ene is None or fg is None:
        return None
    ez = (spec.ez_ene(numbered) if spec.ez_ene else "") or ""
    return (
        f"{ez}{s}-{ene}-{spec.ene_seg[0]}-{fg}-{spec.en_suf}",
        f"{ez}{zs}-{ene}-{spec.ene_seg[1]}-{fg}-{spec.zh_suf}",
    )

@dataclass(frozen=True)
class _Chain:
    kind: str
    en_suf: str            # 词缀后缀: "ol" / "one" / "ene" / "oic acid"
    zh_suf: str            # "醇" / "酮" / "烯" / "酸"
    coda: str              # 饱和词干后接 "an"; 烯/炔词干接 ""
    loc: str | None        # FG 位次字段 (单 locant)
    omit_key: str
    default_omit: bool
    no_loc: str            # "plain"=无位次输出普通名; "none"=返回 None
    omit_rule: object      # (n, loc, omit) -> bool  True=省略
    plain_maps: tuple | None = None     # 俗名表 (en_dict, zh_dict), 查不到时回落生成
    plain_fn: object = None             # 派生命名 (n)->pair — 酸酐从酸派生
    locs: str | None = None             # 多 FG 位次字段 (list) — diol/dione/环二酮
    need: int | None = None             # 多 FG 位次数量要求
    ene_seg: tuple | None = None        # 单烯段 ("en","烯") — 醇/酮
    yne_seg: tuple | None = None        # 炔段 ("yn","炔") — 醇/酮
    mult_seg: dict | None = None        # 多烯段表 {2:("dien","二烯"),...} — 醇
    ene_base: tuple | None = None       # 融合式烯基座 = 单烯后缀 ("enoic acid","烯酸") — 实际后缀 = MULT[m]+基座
    ene_m_min: int = 1                  # 融合式烯最小烯数 (polyene/cyclopolyene 为 2)
    ene_special: object = None          # 融合式单烯俗名钩子 (n, numbered)->pair — amide 丙烯酰胺
    yne_suf: tuple | None = None        # 融合式炔后缀 ("ynoic acid","炔酸") — 酸
    ez_ene: object = None               # (numbered)->str  单烯 E/Z
    ez_ene_multi: object = None         # (numbered)->str  多烯 E/Z
    ene_n_min: int = 4                  # 多烯融合式 n 下限 (polyene 类为 0)
    ene_single_min: int = 2             # 融合式单烯 n 下限 (diacid 为 3)
    ene_omit_aware: bool = False        # 烯段受 omit_ene_locant 影响 (环系 FG)
    cyclic: bool = False                # 恒加环前缀 (cycloalkane/cycloalkene)
    cyclic_unsat: bool = False          # 仅烯/炔段时加环 (cyclopolyene: 无烯回落纯烷烃)
    zh_full: bool = False               # 中文词干保留完整烷烃后缀 "烷" (环烷/回落)
    wrap: object = None                 # (pair, numbered)->pair  整体包裹 (E/Z)
    unsat_polyol: bool = False          # 多 FG 词干支持烯/炔插入 (diol/triol: but-2-ene-1,4-diol)

def _chain_zh_base(n: int) -> str | None:
    z = alkane_zh(n)
    return zh_stem(z) if z else None

def _chain_plain(spec: _Chain, s: str, zs: str, n: int) -> tuple[str, str]:
    """普通名: 俗名表 → 派生命名 → 词干拼接."""
    pair = _pair(*spec.plain_maps, n) if spec.plain_maps else None
    if pair is None and spec.plain_fn is not None:
        pair = spec.plain_fn(n)
    return pair if pair is not None else (f"{s}{spec.coda}{spec.en_suf}", f"{zs}{spec.zh_suf}")

def _chain_polyol_stem(n: int, numbered: dict) -> tuple[str, str] | None:
    """多 FG 不饱和词干: but-2-ene / 丁烷-2-烯 (对齐 _unsat_polyol_names 输出)."""
    enes = numbered.get("ene_locants")
    ene = numbered.get("ene_locant")
    yne = numbered.get("yne_locant")
    if not (enes or ene or yne):
        return None
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    if enes and len(enes) >= 2:
        loc = _pair_loc_str(enes)
        me = {2: "diene", 3: "triene"}.get(len(enes), "")
        mz = {2: "二烯", 3: "三烯"}.get(len(enes), "")
        if not me:
            return None
        return f"{en[:-3]}a-{loc}-{me}", f"{zh}-{loc}-{mz}"
    if yne is not None:
        if numbered.get("omit_yne_locant", False):
            return f"{en[:-3]}yne", f"{zh}炔"
        return f"{en[:-3]}-{yne}-yne", f"{zh}-{yne}-炔"
    if ene is not None:
        if numbered.get("omit_ene_locant", False):
            return f"{en[:-3]}ene", f"{zh}烯"
        return f"{en[:-3]}-{ene}-ene", f"{zh}-{ene}-烯"
    return None

def _chain_names(spec: _Chain, n: int, numbered: dict) -> tuple[str, str] | None:
    """单链词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环; 烯/炔段插入由 spec 数据驱动."""
    alk = _alkane_names(n)
    if not alk:
        return None
    top = _chain_unsat(spec, n, numbered)
    if top is not None:
        if spec.cyclic or spec.cyclic_unsat:
            top = (f"cyclo{top[0]}", f"环{top[1]}")
        return spec.wrap(top, numbered) if spec.wrap is not None else top
    s = _en_stem(n)
    if spec.zh_full:
        zf = alkane_zh(n)
        zs = zf if zf else _chain_zh_base(n)
    else:
        zs = _chain_zh_base(n)
    if s is None or zs is None:
        return None
    if spec.locs is not None:      # 多 FG 位次 (diol/dione/环二酮...): 无位次 → None
        locs = numbered.get(spec.locs)
        if not locs or (spec.need is not None and len(locs) != spec.need):
            return None
        loc_s = ",".join(str(x) for x in locs)
        if spec.unsat_polyol:
            poly_stem = _chain_polyol_stem(n, numbered)
            if poly_stem is not None:
                es, zs2 = poly_stem
                pair = (f"{es}-{loc_s}-{spec.en_suf}", f"{zs2}-{loc_s}-{spec.zh_suf}")
            else:
                pair = (f"{s}{spec.coda}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
        else:
            pair = (f"{s}{spec.coda}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
    else:
        loc = numbered.get(spec.loc) if spec.loc else None
        omit = numbered.get(spec.omit_key, spec.default_omit) if spec.loc else False
        if loc is None:
            if spec.no_loc == "none":
                return None
            pair = _chain_plain(spec, s, zs, n)
        elif spec.omit_rule(n, loc, omit):
            pair = _chain_plain(spec, s, zs, n)
        else:
            pair = (f"{s}{spec.coda}-{loc}-{spec.en_suf}", f"{zs}-{loc}-{spec.zh_suf}")
    if spec.cyclic:
        pair = (f"cyclo{pair[0]}", f"环{pair[1]}")
    return spec.wrap(pair, numbered) if spec.wrap is not None else pair

_AMINE_SPEC = _Chain(kind="amine", en_suf="amine", zh_suf="胺", coda="an",
                     loc="amine_locant", omit_key="omit_amine_locant", default_omit=False,
                     no_loc="plain", omit_rule=_omit_term_locant)

_KIND_TABLE = {
    "alcohol": _Chain(kind="alcohol", en_suf="ol", zh_suf="醇", coda="an",
                      loc="oh_locant", omit_key="omit_oh_locant", default_omit=False,
                      no_loc="plain", omit_rule=_omit_term_locant,
                      plain_maps=(ALCOHOL_EN, ALCOHOL_ZH),
                      ene_seg=("en", "烯"), yne_seg=("yn", "炔"),
                      mult_seg={2: ("dien", "二烯"), 3: ("trien", "三烯"),
                                4: ("tetraen", "四烯"), 5: ("pentaen", "五烯")},
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent),
    "diol": _Chain(kind="diol", en_suf="diol", zh_suf="二醇", coda="ane",
                   loc=None, omit_key="", default_omit=False,
                   no_loc="none", omit_rule=lambda n, loc, omit: False,
                   locs="oh_locants", need=2, zh_full=True, unsat_polyol=True),
    "triol": _Chain(kind="triol", en_suf="triol", zh_suf="三醇", coda="ane",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="none", omit_rule=lambda n, loc, omit: False,
                    locs="oh_locants", need=3, zh_full=True, unsat_polyol=True),
    "ketone": _Chain(kind="ketone", en_suf="one", zh_suf="酮", coda="an",
                     loc="ketone_locant", omit_key="omit_ketone_locant", default_omit=False,
                     no_loc="none",
                     omit_rule=lambda n, loc, omit: n <= 2 and loc == 1,
                     ene_seg=("en", "烯"), yne_seg=("yn", "炔"),
                     ez_ene=ez_for_parent),
    "alkene": _Chain(kind="alkene", en_suf="ene", zh_suf="烯", coda="",
                     loc="ene_locant", omit_key="omit_ene_locant", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: omit or n <= 3,
                     wrap=_with_ez),
    "acid": _Chain(kind="acid", en_suf="oic acid", zh_suf="酸", coda="an",
                   loc=None, omit_key="", default_omit=False,
                   no_loc="plain", omit_rule=lambda n, loc, omit: False,
                   plain_maps=(ACID_EN, ACID_ZH),
                   ene_base=("enoic acid", "烯酸"),
                   yne_suf=("ynoic acid", "炔酸"),
                   ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent),
    "ester": _Chain(kind="ester", en_suf="oate", zh_suf="酸", coda="an",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="plain", omit_rule=lambda n, loc, omit: False,
                    plain_maps=(ESTER_ACYL_EN, ACID_ZH),
                    ene_base=("enoate", "烯酸"),
                    yne_suf=("ynoate", "炔酸"),
                    ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent),
    "thiol": _Chain(kind="thiol", en_suf="thiol", zh_suf="硫醇", coda="ane",
                    loc="sh_locant", omit_key="omit_sh_locant", default_omit=False,
                    no_loc="plain", omit_rule=_omit_term_locant),
    "amine": _AMINE_SPEC,
    "sec_amine": _AMINE_SPEC,
    "tert_amine": _AMINE_SPEC,
    "aldehyde": _Chain(kind="aldehyde", en_suf="anal", zh_suf="醛", coda="an",
                       loc=None, omit_key="", default_omit=False,
                       no_loc="plain", omit_rule=lambda n, loc, omit: False,
                       plain_maps=(ALDEHYDE_EN, ALDEHYDE_ZH),
                       ene_base=("enal", "烯醛"),
                       yne_suf=("ynal", "炔醛"),
                       ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent),
    "nitrile": _Chain(kind="nitrile", en_suf="anenitrile", zh_suf="腈", coda="an",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="plain", omit_rule=lambda n, loc, omit: False,
                      plain_maps=(NITRILE_EN, NITRILE_ZH),
                      ene_base=("enenitrile", "烯腈"),
                      yne_suf=("ynenitrile", "炔腈"),
                      ez_ene=ez_for_parent),
    "polyene": _Chain(kind="polyene", en_suf="diene", zh_suf="二烯", coda="",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="none", omit_rule=lambda n, loc, omit: False,
                      ene_base=("ene", "烯"), ene_m_min=2,
                      ene_n_min=0, wrap=_with_ez),
    "cyclopolyene": _Chain(kind="cyclopolyene", en_suf="", zh_suf="", coda="ane",
                           loc=None, omit_key="", default_omit=False,
                           no_loc="plain", omit_rule=lambda n, loc, omit: False,
                           ene_base=("ene", "烯"), ene_m_min=2,
                           ene_n_min=0, cyclic_unsat=True, zh_full=True),
    "alkyne": _Chain(kind="alkyne", en_suf="yne", zh_suf="炔", coda="",
                     loc="yne_locant", omit_key="omit_yne_locant", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: n in (2, 3),
                     plain_maps=({2: "acetylene", 3: "propyne"}, {2: "乙炔", 3: "丙炔"})),
    "cycloalkane": _Chain(kind="cycloalkane", en_suf="", zh_suf="", coda="ane",
                          loc=None, omit_key="", default_omit=False,
                          no_loc="plain", omit_rule=lambda n, loc, omit: False,
                          cyclic=True, zh_full=True),
    "cycloalkene": _Chain(kind="cycloalkene", en_suf="ene", zh_suf="烯", coda="",
                          loc=None, omit_key="", default_omit=False,
                          no_loc="plain", omit_rule=lambda n, loc, omit: False,
                          cyclic=True),
    "dione": _Chain(kind="dione", en_suf="dione", zh_suf="二酮", coda="ane",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="none", omit_rule=lambda n, loc, omit: False,
                    locs="ketone_locants", need=2),
    "cycloalkanedione": _Chain(kind="cycloalkanedione", en_suf="dione", zh_suf="二酮", coda="ane",
                               loc=None, omit_key="", default_omit=False,
                               no_loc="none", omit_rule=lambda n, loc, omit: False,
                               locs="ketone_locants", need=2, cyclic=True, zh_full=True),
    "cycloalcohol": _Chain(kind="cycloalcohol", en_suf="ol", zh_suf="醇", coda="an",
                           loc="oh_locant", omit_key="omit_oh_locant", default_omit=True,
                           no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                           ene_seg=("en", "烯"), ene_omit_aware=True, cyclic=True),
    "cycloketone": _Chain(kind="cycloketone", en_suf="one", zh_suf="酮", coda="an",
                          loc="ketone_locant", omit_key="omit_ketone_locant", default_omit=True,
                          no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                          ene_seg=("en", "烯"), ene_omit_aware=True, cyclic=True),
    "cycloamine": _Chain(kind="cycloamine", en_suf="amine", zh_suf="胺", coda="an",
                         loc="amine_locant", omit_key="omit_amine_locant", default_omit=True,
                         no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                         cyclic=True),
    "diacid": _Chain(kind="diacid", en_suf="dioic acid", zh_suf="二酸", coda="ane",
                     loc=None, omit_key="", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: False,
                     plain_maps=({2: "oxalic acid"}, {2: "草酸"}),
                     ene_base=("enedioic acid", "烯二酸"),
                     ene_single_min=3,
                     ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent),
    "diamine": _Chain(kind="diamine", en_suf="diamine", zh_suf="二胺", coda="ane",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="none", omit_rule=lambda n, loc, omit: False,
                      locs="amine_locants", need=2, zh_full=True),
    "triamine": _Chain(kind="triamine", en_suf="triamine", zh_suf="三胺", coda="ane",
                       loc=None, omit_key="", default_omit=False,
                       no_loc="none", omit_rule=lambda n, loc, omit: False,
                       locs="amine_locants", need=3, zh_full=True),
    "tetraamine": _Chain(kind="tetraamine", en_suf="tetraamine", zh_suf="四胺", coda="ane",
                         loc=None, omit_key="", default_omit=False,
                         no_loc="none", omit_rule=lambda n, loc, omit: False,
                         locs="amine_locants", need=4, zh_full=True),
    "anhydride": _Chain(kind="anhydride", en_suf="", zh_suf="", coda="ane",
                        loc=None, omit_key="", default_omit=False,
                        no_loc="plain", omit_rule=lambda n, loc, omit: False,
                        plain_fn=_anhydride_from_acid),
    "cycloalkanediol": _Chain(kind="cycloalkanediol", en_suf="diol", zh_suf="二醇", coda="ane",
                              loc=None, omit_key="", default_omit=False,
                              no_loc="none", omit_rule=lambda n, loc, omit: False,
                              locs="oh_locants", need=2, cyclic=True, zh_full=True),
    "amide": _Chain(kind="amide", en_suf="amide", zh_suf="酰胺", coda="an",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="plain", omit_rule=lambda n, loc, omit: False,
                    plain_maps=(AMIDE_EN, AMIDE_ZH),
                    ene_base=("enamide", "烯酰胺"),
                    ene_special=alkenamide_names,
                    yne_suf=("ynamide", "炔酰胺"),
                    ez_ene=_ez_prefix),
}

def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Chain-engine dispatch for table kinds, then specific workers."""
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
        return _chain_names(entry, n, numbered)
  
    if kind == "sulfide":
        return _sulfide_names(numbered)
    if kind == "benzenediol":
        return benzenediol_names(numbered.get("oh_locants"))
    if kind == "phenyl":
        return phenyl_parent_names(numbered)
    if kind == "benzene":
        return benzene_parent_names(numbered)
    if kind == "benzoate":
        return ("benzoate", "苯甲酸")
    top = pyridine_kind_names(kind, numbered, _build_prefix)
    if top is not None:
        return top
    stem = _parent_stem_names(numbered)
    return stem if stem is not None else _alkane_names(n)

def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None

def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)
def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})
from namepredict.layer5.assembler_prefixes import _build_prefix, _prefix_for

def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    from namepredict.layer5.stereo_rs import apply_rs_prefix
    kind, n = _parent_n(numbered)
    effective_kind = _typed_expression_kind(kind, numbered)
    names = _names_for(effective_kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    en, zh = join_kind_name(effective_kind, _prefix_for(numbered, effective_kind, n), names, numbered)
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
