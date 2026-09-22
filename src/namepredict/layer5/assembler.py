"""L5 名称组装主入口：取名后拼前缀、立体与盐后缀。"""
from __future__ import annotations
from dataclasses import replace
import re
from rdkit import Chem
from rdkit.Chem import BondType
from namepredict.constants import (
    ALKOXY_YLOXY_EN, ALKOXY_YLOXY_ZH, AMIDO_RETAINED, AZANE_PAREN_SUF, BIS_EN, BRIDGE_YL_SUFFIX,
    BRIDGE_ZH_YL_SUFFIX, C, ESTER_O_SIDE_KINDS, N, O, OXO_CENTER_KINDS, S,
    BRIDGE_FUSION_YL, CATION_STEMS, MONONUCLEAR_BRIDGE, MONONUCLEAR_HYDRIDES,
    MONONUCLEAR_YL, MONONUCLEAR_ZERO_YL, MULT_EN,
    MULT_ZH, PHOSPHORYL_STEMS, ZH_DIGITS, zh_bridge_root,
)
from namepredict.tools.anchored_table import carbamoyl_prefix_name
from namepredict.tools.re import SUB_LOCANT_RE, alpha_order_key
from namepredict.layer5.chain_engine import (
    _ACYL_HALIDE_BY_HAL, _BENZENE_RETAINED, _KIND_TABLE, _benzene_retained, _chain_names,
)
from namepredict.layer5.stems import (
    _metal_en_prefix, _metal_zh_suffix, join_anion_names,
)
from namepredict.layer5.assembler_prefixes import (
    _ACYL_FRONT_RE, _SIMPLE_CHAIN_YL_RE, _STEREO_LEAD_ENCLOSE_RE, _enclose, _mult_rows,
    _prefix_for, oxo_arm_fence,
)
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

    pairs: list[tuple[str, str]] = []
    for s in subs:
        en, zh = (s.get("en") or "").strip(), (s.get("zh") or "").strip()
        if not en or not zh:
            return None
        oxido = _oxido_arm(s, mol)
        if oxido is not None:  # 酸式 H 已被夺去的 O⁻ 臂：hydroxy → oxido
            en, zh = oxido
        pairs.append((en, zh))
    rows = _mult_rows(pairs, lambda p: p[0], lambda p: p[1], alpha_order_key)
    if len(rows) == 1 and rows[0][2] > 1:  # 同基倍增（P-16.3.2 简单基用 di-）
        en, zh, m = rows[0][:3]
        m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
        if not m_en or not m_zh:
            return None
        return f"{m_en}{en}{stem_en}", f"{m_zh}{zh}{stem_zh}基"
    compound = any("(" in r[0] or "[" in r[0] for r in rows)  # 含自身带括号的复合组分：逐组分连字符 + 方括号围栏（P-16.5.2）
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for i, (en, zh, m, _) in enumerate(rows):
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

_SPIRO_KINDS = ("mono_spiro", "fused_bridged_spiro")  # P-24 螺环 scaffold 直取的 kind（不在 _KIND_TABLE）



def _azane_front_needs_paren(a: dict) -> bool:
    """azane 前端是否须括起再缀 amino（P-63.2.2.1.1）：前端自带多个位次段（复合取代基）。

    单取代前端（2-sulfanylethyl、2-hydroxyethyl）、酰基前端与苯基前端（走 anilino 保留式）
    直接与 amino 融合（…oylamino / …anilino）。
    """
    name = a.get("en") or ""
    if not name or _ACYL_FRONT_RE.search(name) or any(c in name for c in "()[]"):
        return False  # 酰基前端走融合式；已自带括号的前端由 L5 统一升级围栏
    if name.endswith("phenyl"):
        return False  # 苯基前端走 anilino 保留式（P-62.2.1.1），不再单独括起
    return len(SUB_LOCANT_RE.findall(name)) >= 2


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


# free 母体名 → P-29 -yl 取代基形式的双语转换。
_MONONUCLEAR_NAMES = tuple(MONONUCLEAR_HYDRIDES)  # 本模块转换的单核母体氢化物（P-15.4.1 表 2.1）+ 阳离子词干（P-73.1.1）


def _mononuclear_en(en: str) -> str | None:
    """单核氢化物 free 名 → 去氢取代基名。"""
    for en_suf in _MONONUCLEAR_NAMES:
        if en.endswith("-" + en_suf):
            base = en[: -len(en_suf) - 1]
            yl = MONONUCLEAR_YL[en_suf][1]
            if yl in ("sulfonyl", "sulfinyl") and base == "phenyl":  # P-66.4.1 苯+S 桥取保留名 benzenesulfonyl
                return "benzene" + yl
            return base[: -len("phenyl")] + "anilino" if yl == "amino" and base.endswith("phenyl") else base + yl
    return None


def _mononuclear_zh(zh: str) -> str | None:
    """中文组装名去氢：乙基-氧化烷 → 乙氧基等。"""
    for en_suf in _MONONUCLEAR_NAMES:
        zh_suf, _, zy = MONONUCLEAR_YL[en_suf]
        if zh.endswith("-" + zh_suf):
            base = zh[: -len(zh_suf) - 1]
            if zy == "氨基" and base.endswith("苯基"):
                return base[: -len("苯基")] + "苯胺基"
            if en_suf in CATION_STEMS:  # P-73.1.1 阳离子前缀保留烃基尾「基」（甲基铵基，非甲铵基）
                return base + zy
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


def _fused_bridge_name(stem_en: str, a: dict) -> tuple[str, str] | None:
    """双原子桥合一保留前缀：R-亚氨基/硫基 → R-diazenyl/disulfanyl（P-68.3.1.3/.4）。"""
    for (stem, tail_en), (tails_zh, en_suf, zh_suf) in BRIDGE_FUSION_YL.items():
        if stem != stem_en or not a["en"].endswith(tail_en):
            continue
        tail_zh = next((t for t in tails_zh if a["zh"].endswith(t)), None)
        if tail_zh is None:
            continue
        base_zh = a["zh"][: -len(tail_zh)]
        base_en = a["en"][: -len(tail_en)]
        # 前端中文名尾已带「基」（如 甲基二硫代基 去掉尾后为 甲基）：不再补「基」，避免 甲基基…
        zh_head = base_zh if base_zh.endswith("基") else (base_zh + "基" if base_zh else "")
        return base_en + en_suf, zh_head + zh_suf
    return None


def _same_name_groups(subs: list) -> list[tuple[int, str, str]]:
    """相邻同名取代基按数量收拢为 (数量, en, zh)，保持字母序。"""
    out: list[tuple[int, str, str]] = []
    for s in subs:
        if out and out[-1][1] == s["en"]:
            n, en, zh = out[-1]
            out[-1] = (n + 1, en, zh)
        else:
            out.append((1, s["en"], s["zh"]))
    return out


def _thioxo_fused_stem(stem_en: str, stem_zh: str, subs: list[dict]) -> tuple[str, str, list[dict]]:
    """P 锚点无 =O 而恰有一个 =S 时并入词干：<取代基>phosphinothioyl（P-67.1.4.1.1.4）。"""
    if stem_en != "phosphanyl":
        return stem_en, stem_zh, subs
    hit = [s for s in subs if (s.get("en") or "") == "sulfanylidene"]
    if len(hit) != 1:
        return stem_en, stem_zh, subs
    return "phosphinothioyl", "硫代磷酰基", [s for s in subs if s is not hit[0]]


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
        stem_en, stem_zh, subs = _thioxo_fused_stem(stem_en, stem_zh, subs)
        return _phosphoryl_sub_names(subs, stem_en, stem_zh, parent.get("mol"))
    if len(subs) == 1:
        a = subs[0]
        if stem_en == "sulfonyl" and a["en"].endswith("amino") and a["zh"].endswith("氨基"):  # P-66.1.1.4.2：N-取代基与 sulfamoyl 融合
            return a["en"][: -len("amino")] + "sulfamoyl", a["zh"] + "磺酰基"
        if stem_en == "sulfonyl" and a["en"].endswith("anilino") and a["zh"].endswith("苯胺基"):  # P-66.1.1.4.2：N-芳基按苯基并入 sulfamoyl
            bare = a["en"] == "anilino"
            aryl_en = a["en"][: -len("anilino")] + "phenyl"
            aryl_zh = a["zh"][: -len("苯胺基")] + "苯基"
            return ((aryl_en if bare else _enclose(aryl_en)) + "sulfamoyl",
                    _enclose(aryl_zh) + "氨基磺酰基")
        fused = _fused_bridge_name(stem_en, a)  # -N=N-R / -S-S-R：中心与前端合一为 diazenyl / disulfanyl
        if fused is not None:
            return fused
        bridge = _bridge_enclosed_names(a, stem_en, stem_zh)  # 围栏在 L3 一次定形，中英文同步产出，L5 前缀渲染不再二次拆分
        if bridge is not None:
            if (stem_en, stem_zh) != ("azane", "氮烷"):  # N 桥复合前缀仍须 L5 整体围栏（P-16.5.2）
                numbered["bridge_self_enclosed"] = True
            return bridge
        if stem_en == "azane":
            amido = AMIDO_RETAINED.get(a.get("en") or "")  # P-66.1.1.4.3 方法 1：单 N-酰基残基收成 amido 保留式
            if amido is not None:
                return amido
            if _azane_acyl_stereo_lead(a.get("en") or "") or _azane_sub_needs_paren(a):
                w_en = _enclose(a["en"])  # 内层已含括号（立体描述符）时升级为方括号（P-16.5.2 嵌套）
                w_zh = _enclose(a["zh"])
                return f"{w_en}amino", f"{w_zh}氨基"
            if _azane_front_needs_paren(a):
                return f"({a['en']})amino", f"({a['zh']})氨基"  # 复合前端：前端括起再缀 amino（P-63.2.2.1.1）
        en, zh = free_to_yl(f"{a['en']}-{stem_en}", f"{a['zh']}-{stem_zh}", 1,
                            paren=bool(a.get("paren")))[:2]
        return _retained_alkoxy(en, zh) if stem_en == "oxidane" else (en, zh)
    zero = MONONUCLEAR_ZERO_YL.get((stem_en, stem_zh))
    cation = stem_en in CATION_STEMS  # P-73.1.1 阳离子词干：2/3 个烃基臂照样成立（二甲基铵基/三甲基铵基）
    if (not cation and (stem_en, stem_zh) != ("azane", "氮烷")) or zero is None or not 2 <= len(subs) <= (3 if cation else 2):
        return None  # 多取代基仅 N（azane）双烷基/阳离子 N 三烷基成立：O/S 双烷基非标准自由基，明确失败。
    ordered = sorted(subs, key=lambda s: s["en"])
    if len({s["en"] for s in ordered}) == 1:
        base = ordered[0]
        base_zh = base["zh"] if cation else zh_bridge_root(base["zh"])  # 阳离子前缀保留烃基尾「基」（甲基铵基）
        return (f"{MULT_EN[len(ordered)]}{base['en']}{zero[0]}",
                f"{MULT_ZH[len(ordered)]}{base_zh}{zero[1]}")
    aryl = [s for s in ordered if s["en"].endswith("phenyl") and s["zh"].endswith("苯基")]
    if len(aryl) == 1 and not cation:  # P-62.2.1.1：N-芳基-N-某基胺取 anilino（阳离子无此保留名）
        ring = aryl[0]
        other = next(s for s in ordered if s is not ring)
        ring_en, ring_zh = ring["en"][: -len("phenyl")], ring["zh"][: -len("苯基")]
        if alpha_order_key(ring_en) <= alpha_order_key(other["en"]):
            return (f"{ring_en}-N-{other['en']}anilino" if ring_en else f"N-{other['en']}anilino",
                    f"{ring_zh}-N-{other['zh']}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
        return (f"N-{other['en']}-{ring_en}anilino" if ring_en else f"N-{other['en']}anilino",
                f"N-{other['zh']}-{ring_zh}苯胺基" if ring_zh else f"N-{other['zh']}苯胺基")
    groups = _mult_rows(ordered, lambda s: s["en"], lambda s: s["zh"])  # 双不同 N-取代基：首基平铺、其余各基加括号（P-62.2.2.1）；同名基按数量词收拢（P-16.5.1.3.1）
    en_tail = groups[0][0] if groups[0][2] == 1 else MULT_EN[groups[0][2]] + groups[0][0]
    zh_tail = groups[0][1] if groups[0][2] == 1 else MULT_ZH[groups[0][2]] + groups[0][1]
    return (en_tail + "".join(f"{MULT_EN[r[2]]}({r[0]})" if r[2] > 1 else f"({r[0]})" for r in groups[1:]) + zero[0],
            zh_tail + "".join(f"{MULT_ZH[r[2]]}({r[1]})" if r[2] > 1 else f"({r[1]})" for r in groups[1:]) + zero[1])


_STEM_H_PREFIX_RE = re.compile(r"^(?:\d+[a-z]?H[-,])+")  # 词干自带指示氢前缀（1H- / 3H,4H- / 7H-…）


def _split_stem_h_prefix(name: str) -> tuple[str, str]:
    """拆出词干自带的指示氢前缀，返回 (前缀, 余下词干)。"""
    m = _STEM_H_PREFIX_RE.match(name or "")
    return (m.group(0), name[m.end():]) if m else ("", name or "")


def _indicated_h_prefix(parent: dict) -> str:
    """把 L4 定好的指示氢位次拼成前缀（P-58.2.1）。"""
    locants = [str(x) for x in (parent.get("indicated_h_locants") or ()) if str(x)]
    seen: list[str] = []
    for loc in locants:  # 同位次去重：嘌呤 7H/9H 等效位不得双写
        if loc not in seen:
            seen.append(loc)
    return f"{','.join(f'{l}H' for l in seen)}-" if seen else ""


def _with_indicated_h(name: str, ind: str) -> str:
    """把 L4 指示氢前缀落到词干前：词干自带前缀被取代，避免 1H-7H- 双写。"""
    if not ind or name.startswith(ind):
        return name
    return ind + _split_stem_h_prefix(name)[1]


def _h_prefix_atoms(parent: dict, prefix: str) -> list[int]:
    """词干指示氢前缀的位次 → 母体链原子（位次不可映射时返回空表）。"""
    labels = (parent.get("numbering_scaffold") or {}).get("labels")
    chain = list(parent.get("chain") or ())
    if not chain:
        return []
    out: list[int] = []
    for loc in re.findall(r"(\d+[a-z]?)H", prefix or ""):
        for idx, atom in enumerate(chain):  # 缺 labels 时回退链序号（同 atom_locant）
            lab = str(labels[idx]) if labels and len(labels) == len(chain) else str(idx + 1)
            if lab == loc:
                out.append(atom)
                break
    return out


_DIHYDRO_STEM_RE = re.compile(r"^(\d+(?:,\d+)*)-dihydro(.+)$")
_DIHYDRO_STEM_ZH_RE = re.compile(r"^(\d+(?:,\d+)*)-二氢(.+)$")


def _dihydro_stem_kept(parent: dict, locants: str) -> list[str]:
    """滤出词干 'N,M-二氢' 中真正饱和的位次：已带多重键者改由指示氢表达。"""
    mol = parent.get("mol")
    if mol is None:
        return locants.split(",")
    keep: list[str] = []
    for loc in locants.split(","):
        atoms = _h_prefix_atoms(parent, f"{loc}H-")
        if atoms and any(b.GetBondType() != Chem.BondType.SINGLE
                         for b in mol.GetAtomWithIdx(atoms[0]).GetBonds()):
            continue  # 该位已是酮/烯碳，不能再算二氢（P-58.2.1）
        keep.append(loc)
    return keep


def _dihydro_stem_fix(parent: dict, stem_en: str, stem_zh: str) -> tuple[str, str]:
    """词干自带的 'N,M-二氢' 中已带多重键的位次改为指示氢（P-58.2.1）：
    5-oxo-2,5-dihydrofuran-3-yl → 5-oxo-2H-furan-3-yl（同一位不得既是酮又是 CH2）。"""
    m = _DIHYDRO_STEM_RE.match(stem_en or "")
    if not m:
        return stem_en, stem_zh
    locs = m.group(1)
    keep = _dihydro_stem_kept(parent, locs)
    if len(keep) == len(locs.split(",")):
        return stem_en, stem_zh
    rest_en, rest_zh = m.group(2), stem_zh
    mz = _DIHYDRO_STEM_ZH_RE.match(stem_zh or "")
    if mz:
        rest_zh = mz.group(2)
    if not keep:  # 全部位次都带多重键：只剩母体名，由 L4 指示氢另行补位
        return rest_en, rest_zh
    return f"{','.join(keep)}H-{rest_en}", f"{','.join(keep)}H-{rest_zh}"


def _nh_only(parent: dict, ind: str) -> bool:
    """指示氢是否全落在环内氮：仅 NH 才在非 forced 路径补前缀（P-58.2.1）。"""
    mol = parent.get("mol")
    atoms = _h_prefix_atoms(parent, ind)
    return bool(atoms) and mol is not None and all(
        mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in atoms if i < mol.GetNumAtoms())


def _ring_n_free(parent: dict) -> bool:
    """母体环内氮是否都未被取代：N-取代基环（如 N-糖基尿嘧啶）不再写 NH 前缀。"""
    mol = parent.get("mol")
    chain = list(parent.get("chain") or ())
    if mol is None or not chain:
        return False
    ring = set(chain)
    for i in chain:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() != 7:
            continue
        if any(n.GetIdx() not in ring for n in atom.GetNeighbors()):
            return False
    return True


def _stem_prefix_stale(parent: dict, prefix: str) -> bool:
    """词干自带指示氢位次是否已成为羰基碳（P-58.2.1）：如茚-1,3-二酮的 C1。"""
    atoms = _h_prefix_atoms(parent, prefix)
    ring = set(parent.get("chain") or ())
    mol = parent.get("mol")
    if not atoms or mol is None:
        return False
    for i in atoms:
        if i >= mol.GetNumAtoms():
            return False
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() > 0:
            return False
        if not any(b.GetBondType() != Chem.BondType.SINGLE
                   and b.GetOtherAtomIdx(i) not in ring
                   and mol.GetAtomWithIdx(b.GetOtherAtomIdx(i)).GetAtomicNum() != 6
                   for b in atom.GetBonds()):
            return False
    return True


def _ensure_parent_stem(numbered: dict) -> bool:
    """母体词干注入：螺环 P-24 → 桥环 von Baeyer → 稠环稠合 → 大环杂单环生成式。"""
    parent = numbered.get("parent") or {}
    if parent.get("stem_en") and parent.get("stem_zh"):
        return True  # 已注入则幂等返回；各子注入器共用此守卫
    if parent.get("fbs_node") is not None:
        return _ensure_fbs_stem(numbered)
    if parent.get("spiro_node") is not None:
        return _ensure_spiro_stem(numbered)
    if parent.get("bridged_node") is not None:
        return _ensure_bridged_stem(numbered)
    if not _ensure_fused_stem(numbered):
        return False
    return _ensure_generated_stem(numbered)


def _ensure_generated_stem(numbered: dict) -> bool:
    """P-23.3.1 生成式词干：>10 元纯杂单环无保留名"""
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), list(parent.get("chain") or ())
    if mol is None or parent.get("scaffold_id") != "carbocycle" or len(chain) <= 10:
        return True  # ≤10 元环归 Hantzsch-Widman 保留名
    from namepredict.layer5.skeleton_replacement import prefix_from_chain
    a_en, a_zh = prefix_from_chain(mol, chain)
    if not a_en:
        return a_en is not None  # 无杂原子（纯环烷）放行；词表外元素判失败
    from namepredict.layer5.stems import alkane_en, alkane_zh, stem_forms
    en, zh = alkane_en(len(chain)), alkane_zh(len(chain))
    if not en or not zh:
        return False
    (full_en, full_zh), (bare_en, bare_zh) = stem_forms(f"{a_en}cyclo", f"{a_zh}环", en, zh)
    parent["stem_en"], parent["stem_zh"] = full_en, full_zh
    parent["stem_bare_en"], parent["stem_bare_zh"] = bare_en, bare_zh
    parent["stem_generated"] = True  # 'a' 前缀生成式词干：位次恒显式，自由基取裸词干
    return True


def _ensure_spiro_stem(numbered: dict) -> bool:
    """P-24 螺环词干注入。"""
    from namepredict.layer5.spiro_namer import spiro_parent_names
    return _ensure_ring_stem(numbered, "spiro_node", spiro_parent_names)


def _ensure_bridged_stem(numbered: dict) -> bool:
    """P-23 桥环词干注入。"""
    from namepredict.layer5.bridged_namer import bridged_parent_names
    return _ensure_ring_stem(numbered, "bridged_node", bridged_parent_names)


def _ensure_ring_stem(numbered: dict, node_key: str, namer) -> bool:
    """P-23/P-24 通用环词干注入：饱和取完整名、带不饱和取裸词干（同 _exo_ring_spec 惯例）。"""
    parent = numbered.get("parent") or {}
    mol, node = parent.get("mol"), parent.get(node_key)
    chain = list(parent.get("chain") or ())
    if mol is None or node is None or not chain:
        return False
    names = namer(mol, node, chain)
    if names is None:
        return False
    (full_en, full_zh), (bare_en, bare_zh) = names
    parent["stem_bare_en"], parent["stem_bare_zh"] = bare_en, bare_zh
    unsat = bool(numbered.get("ene_locants") or numbered.get("yne_locants"))
    parent["stem_en"], parent["stem_zh"] = (bare_en, bare_zh) if unsat else (full_en, full_zh)
    return True


def _ensure_fbs_stem(numbered: dict) -> bool:
    """P-24.5~24.7 组分式螺环词干注入：组分名自带不饱和，故不设裸词干形态。"""
    parent = numbered.get("parent") or {}
    mol, node = parent.get("mol"), parent.get("fbs_node")
    if mol is None or node is None:
        return False
    from namepredict.layer5.spiro_namer import fbs_parent_name
    name = fbs_parent_name(mol, node)
    if name is None or not name[0] or not name[1]:
        return False
    parent["stem_en"], parent["stem_zh"] = name
    return True


def _ensure_fused_stem(numbered: dict) -> bool:
    """未注册稠环词干注入：由 fused_tree 组装稠合 base 名。"""
    parent = numbered.get("parent") or {}
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

def _oxoacid_arm_zh(zh: str) -> str:
    """含氧酸 O-侧臂中文词：简单基去「基」，复合名原样保留。"""
    if not zh.endswith("基") or "-" in zh or zh.startswith("("):
        return zh
    stem = zh[:-1]
    if len(stem) >= 2 and all(c in ZH_DIGITS + "十" for c in stem):
        return f"{stem}烷基"
    return stem


def _fenced_arm(name: str, sub: dict, mol) -> str:
    """O-侧臂名围栏：判据命中时整体加括号（内含圆括号则升为方括号）。"""
    return _enclose(name) if oxo_arm_fence(name, sub, mol) else name


def _fenced_arms(arms: list[dict], mol) -> list[dict]:
    """按围栏判据改写 O-侧臂的双语名。"""
    return [{**a, "en": _fenced_arm(a.get("en") or "", a, mol),
             "zh": _fenced_arm(a.get("zh") or "", a, mol)} for a in arms]


def _fenced_arm_phosphate(name: str, sub: dict, mol) -> str:
    """磷酸酯 O-侧臂围栏（P-67.1.3）：环上复合臂或带多个位次的支链臂整体括起。"""
    if not name or "(" in name or name[:1] == "[":
        return name  # 已自带括号的复合名不再二次围栏
    if not sub.get("paren"):
        return name  # 简单基（methyl/ethyl/…）无歧义，平铺
    hit = oxo_arm_fence(name, sub, mol) or sum(
        1 for seg in name.split("-") if seg[:1].isdigit()) > 1
    return _enclose(name) if hit else name


_P_ACYL_ARM_TAIL = ("phosphoryl", "phosphanyl", "phosphinothioyl")  # P 酰基臂名尾（含 P=S 词干）


def _is_condensed_phosphate(arms: list[dict]) -> bool:
    """侧臂中是否含缩合磷酸的 P 酰基臂（-O-P(=O)(OH)- 的 phosphoryl 名），即 P-O-P 二酯。"""
    return any((a.get("en") or "").endswith(_P_ACYL_ARM_TAIL) for a in arms)


def _fenced_arms_phosphate(arms: list[dict], mol) -> list[dict]:
    """按磷酸酯专用判据改写 O-侧臂的双语名。"""
    if _is_condensed_phosphate(arms) and any((a.get("en") or "")[:1] == "(" for a in arms):
        # 缩合磷酸酯（P-O-P）中任一侧臂以前导立体描述符起（非完整围栏）：两臂须同时整体围栏，否则臂界不清（P-16.5.1.3.1/P-67.1.3）
        return [{**a, "en": _enclose(a.get("en") or ""), "zh": _enclose(a.get("zh") or "")}
                for a in arms]
    if _all_arms_stereo_lead(arms):  # 全部臂名均带前导立体描述符：须整体围栏（P-16.5.1.3.1）
        return [{**a, "en": _enclose(a.get("en") or ""), "zh": _enclose(a.get("zh") or "")}
                for a in arms]
    return [{**a, "en": _fenced_arm_phosphate(a.get("en") or "", a, mol),
             "zh": _fenced_arm_phosphate(a.get("zh") or "", a, mol)} for a in arms]


def _all_arms_stereo_lead(arms: list[dict]) -> bool:
    """多个 O-侧臂是否全部以前导立体描述符开头（单臂不受 P-16.5.1.3.1 的第二臂规则支配）。"""
    ens = [a.get("en") or "" for a in arms]
    return len(ens) >= 2 and all(_STEREO_LEAD_ENCLOSE_RE.match(e) for e in ens)


def join_oxoacid_name(pre: tuple[str, str], names: tuple[str, str], numbered: dict) -> tuple[str, str] | None:
    """含氧酸中心母体整名（P-67.1.3）：取代前缀 + O-侧臂 + 词尾 + 金属盐。"""
    tail_en = join_parent_name(pre[0], names[0])
    tail_zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    parent = numbered.get("parent") or {}
    arms_in = _o_side_arms(numbered)
    if parent.get("oxo_kind") == "phosphate":  # 磷酸酯臂按专用判据加围栏（P-67.1.3）
        arms_in = _fenced_arms_phosphate(arms_in, parent.get("mol"))
    else:  # 其余中心母体按判据加围栏
        arms_in = _fenced_arms(arms_in, parent.get("mol"))
    arms = _join_o_side_arms(arms_in, group=True, arm_zh_fn=_oxoacid_arm_zh)
    if arms is None:
        return None
    alk_en, alk_zh = arms
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


def _c1_amino(parent: dict) -> bool:
    """单碳母体（C1 保留名）的官能团碳是否直连 N（氨基甲酸酯类，P-66.3.2）。"""
    mol, chain = parent.get("mol"), list(parent.get("chain") or ())
    if mol is None or len(chain) != 1:
        return False
    atom = mol.GetAtomWithIdx(int(chain[0]))
    return atom.GetAtomicNum() == 6 and any(
        n.GetAtomicNum() == 7 for n in atom.GetNeighbors())


def _c1_two_hetero(parent: dict) -> tuple[list[int], list[int], dict[int, int]] | None:
    """单碳母体官能团碳的 (单键 N, 单键 O, 双键杂原子→索引)；非单碳碳中心返回 None。"""
    mol, chain = parent.get("mol"), list(parent.get("chain") or ())
    if mol is None or len(chain) != 1:
        return None
    atom = mol.GetAtomWithIdx(int(chain[0]))
    if atom.GetAtomicNum() != C:
        return None
    sgl_n: list[int] = []
    sgl_o: list[int] = []
    dbl: dict[int, int] = {}
    for nb in atom.GetNeighbors():
        z = nb.GetAtomicNum()
        if z <= 1:
            continue  # 氢与锚定哑原子不计
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
        if bond.GetBondType() == BondType.SINGLE:
            (sgl_n if z == N else sgl_o if z == O else []).append(nb.GetIdx())
        elif bond.GetBondType() == BondType.DOUBLE:
            dbl[z] = nb.GetIdx()
        else:
            return None
    return sgl_n, sgl_o, dbl


def _is_dbl_hetero_sub(mol, sub: dict, het_idx: int) -> bool:
    """取代基块是否即母体碳上那个双键杂原子（=S/=N 已并入保留名，须撤下）。"""
    return het_idx in (sub.get("atoms") or ())


def _sub_on_n(mol, sub: dict, n_atoms: set[int]) -> int | None:
    """取代基块在母体侧所连的氮原子索引；不连氮返回 None。"""
    for a in sub.get("atoms") or ():
        for nb in mol.GetAtomWithIdx(int(a)).GetNeighbors():
            if nb.GetIdx() in n_atoms:
                return nb.GetIdx()
    return None


def _urea_subs(numbered: dict, n_atoms: list[int], locants: tuple[str, str],
               force: bool = False) -> str | None:
    """脲/硫脲/胍母体：N-取代基按数字位次引用（P-14.4 最低位次给先引用者）。"""
    mol = (numbered.get("parent") or {}).get("mol")
    subs = numbered.get("substituents") or []
    if mol is None:
        return None
    groups: dict[int, list] = {}
    for s in subs:
        n_idx = _sub_on_n(mol, s, set(n_atoms))
        if n_idx is not None:
            groups.setdefault(n_idx, []).append(s)
    for s in subs:
        s["kind"] = "side"  # 改走数字位次通道，不再按 N- 前缀渲染
    if len(subs) < 2 and not force:  # 全分子仅一个取代基、无歧义：位次省略（(4-甲基苯基)脲）
        return None
    order = sorted(groups, key=lambda i: alpha_order_key(groups[i][0].get("en") or ""))
    for pos, n_idx in enumerate(order[: len(locants)]):
        for s in groups[n_idx]:
            s["locant"] = locants[pos]
    return "urea"


def _carbamothioylamino_prefix(side_en: str, side_zh: str) -> tuple[str, str] | None:
    """N-侧胺名 → <R>carbamothioyl / <R>氨基硫代羰基（P-66.1.1.4 硫代氨基甲酸残基；尾「amino/氨基」由自由价渲染补出）。"""
    if not side_en.endswith("amino") or not side_zh.endswith("氨基"):
        return None
    stem_en, stem_zh = side_en[: -len("amino")], side_zh[: -len("氨基")]
    if not stem_en or not stem_zh:
        return None
    if not stem_zh.endswith("基"):  # 中文胺名去「氨基」会连「基」一起削掉（环己氨基 → 环己），须补回
        stem_zh += "基"
    return _enclose(stem_en) + "carbamothioyl", _enclose(stem_zh) + "氨基硫代羰基" + "基"


def _c1_retained(numbered: dict) -> tuple[str, str] | None:
    """单碳母体带两个杂原子时的保留名（P-66.3 脲/硫脲/胍，P-65.2.1.5 carbamoyl）。"""
    parent = numbered.get("parent") or {}
    env = _c1_two_hetero(parent)
    if env is None:
        return None
    sgl_n, sgl_o, dbl = env
    kind = parent.get("kind")
    mol = parent.get("mol")
    if kind == "radical" and len(sgl_n) == 1 and dbl.get(O) is not None and not sgl_o:
        n_idx = sgl_n[0]
        side = next((s for s in (numbered.get("substituents") or []) if n_idx in (s.get("atoms") or ())), None)
        if side is None:  # 伯酰胺：无可引用 N-取代基
            parent["subs_consumed"] = True
            return "carbamoyl", "氨基甲酰基"
        named = carbamoyl_prefix_name(side.get("en") or "", side.get("zh") or "",
                                      in_ring=mol.GetAtomWithIdx(n_idx).IsInRing())
        if named is None:
            return None
        parent["subs_consumed"] = True
        return named
    if kind == "radical" and len(sgl_n) == 1 and dbl.get(S) is not None and not sgl_o:
        n_idx = sgl_n[0]  # 硫代氨基甲酸残基：=S 并入 carbamothioyl，N-侧取代基并入前缀（P-66.1.1.4）
        side = next((s for s in (numbered.get("substituents") or []) if n_idx in (s.get("atoms") or ())), None)
        if side is None:
            return None  # 伯硫代酰胺（-C(=S)NH₂）：无 N-取代基可并入
        named = _carbamothioylamino_prefix(side.get("en") or "", side.get("zh") or "")
        if named is None:
            return None
        parent["subs_consumed"] = True
        return named
    if kind == "amide" and len(sgl_n) == 2 and dbl.get(O) is not None:
        parent["locant_kind"] = _urea_subs(numbered, sgl_n, ("1", "3")) and "urea"
        return "urea", "脲"
    if kind == "amine" and len(sgl_n) == 2 and dbl.get(S) is not None:
        numbered["substituents"] = [s for s in (numbered.get("substituents") or [])
                                    if not _is_dbl_hetero_sub(mol, s, dbl[S])]
        parent["locant_kind"] = _urea_subs(numbered, sgl_n, ("1", "3")) and "thiourea"
        return "thiourea", "硫脲"
    if kind == "amine" and len(sgl_n) == 2 and dbl.get(N) is not None:
        numbered["substituents"] = [s for s in (numbered.get("substituents") or [])
                                    if not _is_dbl_hetero_sub(mol, s, dbl[N])]
        parent["locant_kind"] = _urea_subs(numbered, sgl_n, ("2", "3"), force=True) and "guanidine"
        return "guanidine", "胍"
    if kind == "ester" and len(sgl_o) == 2 and dbl.get(O) is not None and not sgl_n:
        return "carbonate", "碳酸"  # P-65.6.3.3 碳酸二酯：O-侧臂由酯拼接消费
    return None


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """链引擎按表 kind 派发，再转具体 worker。"""
    parent = numbered.get("parent") or {}
    retained = _c1_retained(numbered)  # 单碳双杂原子保留名优先于系统名
    if retained is not None:
        return retained
    if kind == "radical" and parent.get("radical_anchor_element"):
        return _mononuclear_radical_names(numbered)
    entry = _KIND_TABLE.get(kind)
    if kind == "acyl_halide":
        entry = _ACYL_HALIDE_BY_HAL.get(parent.get("hal_z")) or entry
    if parent.get("stem_bare_en") and kind in ("bridged", "alkane", *_SPIRO_KINDS):
        # 无主特征基团的环系（螺环/桥环 / >10 元杂单环）：须走链引擎才会渲染 ene/yne 位次
        return _chain_names(replace(_KIND_TABLE["alkane"],
                                    stem=(parent["stem_bare_en"], parent["stem_bare_zh"]),
                                    coda=""), n, numbered)
    if entry is not None:
        if kind == "ester" and parent.get("thio_side"):  # P-65.6.3.3.7.1 硫代羧酸 S-酯：thioate/硫酯词尾，无 C1/C2 保留名
            entry = replace(entry, coda="ane", en_suf="thioate", zh_suf="硫",
                            ene_base=("enethioate", "烯硫"), yne_suf=("ynethioate", "炔硫"), variant=None)
        sid = parent.get("scaffold_id")
        if kind == "alkane" and parent.get("fused_tree") and sid != "benzene":  # 未注册稠环无 FG：词干注入已完成，返回稠合 base 名
            return _parent_stem_names(numbered)
        if sid == "benzene" and kind == "alkane":  # 苯 base：母体名由 sid 驱动
            return ("benzene", "苯")
        if parent.get("stem_en") and parent.get("stem_zh"):  # 环式 FG 的 locant omit 由 L4 决定，coda 置空
            stem = _dihydro_stem_fix(parent, parent["stem_en"], parent["stem_zh"])
            if kind == "radical" and parent.get("stem_generated"):  # 生成式词干（'a' 前缀大环）的自由基取裸词干（P-29.2）：…tetrazacyclododec-1-yl
                stem = (parent["stem_bare_en"], parent["stem_bare_zh"])
            entry = replace(entry, stem=stem, coda="",
                            omit_rule=lambda n, loc, omit: bool(omit), aromatic=(sid == "benzene"))
        elif sid == "carbocycle":  # 单环饱和烃自由基按 P-29.2 省略 1 位自由价，余沿用 L4 omit
            rule = (lambda n, loc, omit: loc == 1) if kind == "radical" \
                else (lambda n, loc, omit: bool(omit))
            entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=rule)
        # 苯环单 FG 取 scaffold 专属保留名（P-61.2 表）；环外硫代羧酸 S-酯走 carbothioate 词干。
        if sid == "benzene" and kind == "ester" and parent.get("thio_side"):
            sc_variant = _benzene_retained("benzenecarbothioate", "苯硫代甲酸")
        else:
            sc_variant = _BENZENE_RETAINED[kind] if sid == "benzene" and kind in _BENZENE_RETAINED \
                else (entry.variant or {}).get(sid)
        if sc_variant is not None:
            entry = replace(entry, variant=sc_variant)
        result = _chain_names(entry, n, numbered)
        if result and _c1_amino(parent):
            if kind == "ester":  # P-66.3.2 氨基甲酸酯：酸碳连 N → carbamate/氨基甲酸
                result = (result[0].replace("formate", "carbamate"),
                          result[1].replace("甲酸", "氨基甲酸"))
            elif kind == "acid":  # P-66.3.2 氨基甲酸：酸碳连 N → carbamic acid/氨基甲酸
                result = (result[0].replace("formic acid", "carbamic acid"),
                          result[1].replace("甲酸", "氨基甲酸"))
                parent["locant_kind"] = "carbamic_acid"  # P-65.2.1.1 N-取代基不带 N 位次（dimethylcarbamic acid）
                for s in numbered.get("substituents") or []:
                    s["kind"] = "side"  # 撤下 N- 前缀通道，改由位次省略渲染
            elif kind == "acyl_halide":  # P-66.1.1.4.1 氨基甲酰卤：carbamoyl/氨基甲酰
                result = (result[0].replace("formyl", "carbamoyl"),
                          result[1].replace("甲酰", "氨基甲酰"))
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
def join_parent_name(prefix: str, parent: str) -> str:
    """拼接前缀与母体名（数字/1H- 前导时加连字符）。"""
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    needs_hyphen = bool(stem) and (stem[0].isdigit() or stem[0] == "[" or stem.startswith("1H-"))
    body = f"{prefix}-{stem}" if needs_hyphen else f"{prefix}{stem}"
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
        named = [s for s in arms if (s.get("en") or "").strip()]
        rows = _mult_rows(named, lambda s: (s.get("en") or "").strip(),
                          lambda s: arm_zh_fn(s.get("zh") or ""),
                          alpha_order_key)  # P-14.5：O-侧臂同按字母数字序引用
        parts_en: list[str] = []
        parts_zh: list[str] = []
        for en, zh, m, _ in rows:
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


def _fenced_arm_en(name: str) -> str:
    """O 侧臂英文围栏：前导位次或自带括号的复合名须括起。"""
    if not name or not (name[0].isdigit() or "(" in name):
        return name
    return f"[{name}]" if "(" in name else f"({name})"


def _fenced_arm_zh(name: str) -> str:
    """O 侧臂中文围栏：前导位次或自带括号的复合名须括起。"""
    if not name or not (name[0].isdigit() or "(" in name or "[" in name):
        return name
    return f"[{name}]" if ("(" in name or "[" in name) else f"({name})"


def _is_thio_side(numbered) -> bool:
    """母体是否为硫代羧酸 S-酯（P-65.6.3.3.7.1）：S 侧臂改用 thioate 词尾与斜体 S。"""
    return bool((numbered or {}).get("parent", {}).get("thio_side"))


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str] | None:
    """拼接酯名：O 侧作前缀、酸侧作主体（P-16.3.2 倍增）；S 侧走 thioate 整名。"""
    en, zh = names
    thio = _is_thio_side(numbered)
    arms_in = _o_side_arms(numbered)
    if thio:
        arms_in = [{**a, "en": _fenced_arm_en(a.get("en") or "")} for a in arms_in]
    arms = _join_o_side_arms(arms_in, group=False, arm_zh_fn=_fenced_arm_zh if thio else _zh_alkoxy_part)
    # print("_join_o_side_arms",arms,numbered)
    if arms is None:
        return None
    alk_en, alk_zh = arms
    mid = join_parent_name(pre_en, en)
    if thio:  # S-乙基 辛硫酯：位次符号在最前，酯词尾接在母体名后
        return (f"S-{alk_en} {mid}" if alk_en else f"S-{mid}",
                f"S-{alk_zh}{join_parent_name(pre_zh, zh)}酯" if alk_zh else join_parent_name(pre_zh, zh))
    en = f"{alk_en} {mid}" if alk_en else mid
    zh = f"{join_parent_name(pre_zh, zh)}{alk_zh}酯"
    return en, zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str] | None:
    """按 kind 分派：O-侧臂母体（酯/磷酸）走 O-侧拼接，其余走普通母体拼接。"""
    if kind in ESTER_O_SIDE_KINDS:
        if kind in OXO_CENTER_KINDS:  # 中心自任母体（磷酸/膦酸/硫酸酯）：臂 + 功能母体词尾
            return join_oxoacid_name(pre, names, numbered)
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
    ind = _indicated_h_prefix(parent)  # L4 位次为准（P-58.2.1）：词干自带前缀被取代，避免 1H-7H- 双写
    en, zh = names
    if ind and (pre[0] or parent.get("indicated_h_forced") or parent.get("hydro_fallback")
                or (_nh_only(parent, ind) and _ring_n_free(parent))):  # 未取代 NH 环补指示氢
        en, zh = _with_indicated_h(en, ind), _with_indicated_h(zh, ind)
    elif not ind:  # 无 L4 位次时，词干自带前缀若已无 H（如茚二酮 C1）则失效
        s_en, s_zh = _split_stem_h_prefix(en), _split_stem_h_prefix(zh)
        if s_en[0] and _stem_prefix_stale(parent, s_en[0]):
            en, zh = s_en[1], s_zh[1]
    if not pre[0]:
        return en, zh
    return join_parent_name(pre[0], en), join_parent_name(pre[1], zh)


def _cation_insert(en: str, stem: str, base: str) -> tuple[int, int] | None:
    """环阳离子后缀插入点：(切点, 吞字)；跳过连接元音，免切出 iuma 乱词。"""
    at = en.find(stem)
    if at < 0:
        return None
    pos = at + len(stem)
    if stem == base and en[pos:pos + 1] == "a":
        pos += 1  # 连接元音：heptadec → heptadeca-…
    return (pos - 1, 1) if en[pos - 1:pos] == "e" else (pos, 0)


def _zh_ring_cation(zh: str, stem_zh: str, suffix: str) -> str:
    """中文环阳离子：母体词干后插 -{位次}-正离子（中化会 6.7.2，不用“鎓”）。"""
    if not zh or not stem_zh or "鎓" in zh:
        return zh
    for stem in (stem_zh, _split_stem_h_prefix(stem_zh)[1]):  # 词干自带 1H- 时按去前缀词干定位
        if stem and stem in zh:
            at = zh.index(stem) + len(stem)
            return zh[:at] + suffix + zh[at:]
    return zh


def join_ring_cation_suffix(numbered: dict, names: tuple[str, str]) -> tuple[str, str]:
    """环内 N+/O+ → 母体名缀 -{位次}-ium/-{位次}-鎓（P-62.4.1）。"""
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
        base, ium, zh_suf = stem_en[:-1], f"{stem_en[:-3]}enylium", "鎓"
    elif stem_en.endswith("e"):
        base, ium, zh_suf = stem_en[:-1], f"{stem_en[:-1]}-{loc}-ium", f"-{loc}-鎓"
    else:
        base, ium, zh_suf = stem_en, f"{stem_en}-{loc}-ium", f"-{loc}-鎓"
    zh = _zh_ring_cation(zh, parent.get("stem_zh") or "", zh_suf)
    if stem_en.endswith("ene"):  # 色烯型氧鎓保留名：chromene → chromenylium
        token = base + "e" if base + "e" in en else base
        return (en.replace(token, ium, 1), zh) if token in en else names
    hit = _cation_insert(en, stem_en, base) if stem_en in en else _cation_insert(en, base, base)
    if hit is None:
        return names
    cut, drop = hit  # 其余在母体词干后插入 -{位次}-ium
    return en[:cut] + f"-{loc}-ium" + en[cut + drop:], zh

def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    """组装入口：取名 → 前缀 → 阴离子/R-S/金属盐后缀。"""
    from namepredict.layer5.stereo import join_ez_prefix, join_rs_prefix
    # print("layer5 assembling!!!")
    parent = numbered.get("parent") or {}  # 母体 kind 与碳数 n（无则 0）
    kind, n = parent.get("kind"), int(parent.get("n_carbons") or 0)
    if not _ensure_parent_stem(numbered):
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
