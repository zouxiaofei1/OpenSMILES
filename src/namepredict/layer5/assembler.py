"""L5 名称组装主入口：由链引擎取名后按 kind 拼接前缀、立体（E/Z、R/S）与盐类后缀。"""
from __future__ import annotations
from dataclasses import replace

from namepredict.layer5.chain_engine import _KIND_TABLE, _alkane_names, _chain_names
from namepredict.layer5.stems import maybe_anion_names, maybe_metal_salt_names
from namepredict.layer5.assembler_prefixes import _prefix_for
from namepredict.layer5.stereo import _split_stereo_lead as _stereo_lead
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult:
    """构造失败 NameResult（success=False，meta 供诊断）。"""
    return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})


def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult:
    """构造成功 NameResult（success=True，记录耗时与来源）。"""
    return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

# 稠环/杂环词干（en_stem, zh_stem, aromatic）：主路径按 scaffold_id 注入词干；
# aromatic 环上的醇统一用"酚"（苯酚系，含萘/吡啶/吲哚/喹啉），非 variant 特例。
_RING_STEM = {
    "naphthalene": ("naphthalen", "萘", True),
    "indole": ("indol", "吲哚", True),
    "pyridine": ("pyridin", "吡啶", True),
    "quinoline": ("quinolin", "喹啉", True),
    "benzene": ("benzene", "苯", True),
}

# 稠环/杂环完整 base 名（-carboxylic acid 用完整词干，如 naphthalene-1-carboxylic acid）。
_RING_BASE = {
    "naphthalene": ("naphthalene", "萘"),
    "quinoline": ("quinoline", "喹啉"),
    "pyridine": ("pyridine", "吡啶"),
}


def _scaffold_id(numbered: dict) -> str | None:
    """取母体的 scaffold_id（carbocycle/benzene/稠环…）。"""
    return (numbered.get("parent") or {}).get("scaffold_id")


def _exocyclic_acid_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """环酸（carbocycle/稠环 exocyclic COOH）→ cyclohexanecarboxylic acid 式系统名。"""
    parent = numbered.get("parent") or {}
    facts = parent.get("principal_expression_facts")
    sid = parent.get("scaffold_id")
    if not facts or facts.group_class.value != "acid" or facts.relation.value != "exocyclic":
        return None
    if facts.multiplicity != 1:
        return None
    if sid == "carbocycle":
        base = _alkane_names(n)
        return (f"cyclo{base[0]}carboxylic acid", f"环{base[1]}羧酸") if base else None
    if sid in _RING_BASE:
        base = _RING_BASE[sid]
        # 羧基位次：L4 已算出的酸 locant；无则默认省略（1 位）。
        rec = next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == "acid"), None)
        locs = rec.get("locants") if rec else None
        if locs:
            loc = ",".join(str(x) for x in locs)
            return (f"{base[0]}-{loc}-carboxylic acid", f"{base[1]}-{loc}-羧酸")
        return (f"{base[0]}carboxylic acid", f"{base[1]}羧酸")
    return None


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """链引擎按表 kind 派发，再转具体 worker。"""
    if kind == "acid":
        exo = _exocyclic_acid_names(n, numbered)
        if exo:
            return exo
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
        sid = _scaffold_id(numbered)
        if sid == "benzene" and kind == "alkane":
            # 苯 base：无主 FG 的苯，母体名由 sid 驱动（L2 已把纯苯 kind 收敛为 alkane）。
            return ("benzene", "苯")
        if sid in _RING_STEM:
            # 环式 FG 的 locant omit 由 L4 算出的 omit 标志决定；aromatic 标记随 spec 传递，由 _chain_names 消费（芳香醇→酚）。
            stem_en, stem_zh, aromatic = _RING_STEM[sid]
            entry = replace(entry, stem=(stem_en, stem_zh), omit_rule=lambda n, loc, omit: bool(omit),
                            aromatic=aromatic)
        elif sid == "carbocycle":
            entry = replace(entry, cyclic=True, ene_loc_omit=True,
                            omit_rule=lambda n, loc, omit: bool(omit))
        # 苯环单 FG → scaffold 专属保留名 variant (phenol/benzoic…); 开链取 None 键 (acid 草酸)。
        sc_variant = (entry.variant or {}).get(sid)
        if sc_variant is not None:
            entry = replace(entry, variant=sc_variant)
        return _chain_names(entry, n, numbered)

    if kind == "phenyl":
        return "phenyl", "苯基"
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
def join_parent_name(prefix: str, parent: str) -> str:
    """拼接前缀与母体名（数字/1H- 前导时加连字符）。"""
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if stem[:1].isdigit() or stem.startswith("1H-") else f"{prefix}{stem}"
    return f"{stereo}{body}"

def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    """拼接酯名：O 侧烷基作前缀、酸侧作主体（XX 酸YY酯）。"""
    en, zh = names
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    alk_en, alk_zh = o[0].get("en") or "", (o[0].get("zh") or "").rstrip("基")
    st, body = _stereo_lead(en)
    mid = f"{pre_en}{body}" if pre_en else body
    en = f"{alk_en} {st}{mid}" if alk_en else f"{st}{mid}"
    stz, bodyz = _stereo_lead(zh)
    midz = f"{stz}{pre_zh}{bodyz}" if pre_zh else f"{stz}{bodyz}"
    zh = f"{midz}{alk_zh}酯"
    return en, zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str]:
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


def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    """组装入口：取名 → 前缀 → 阴离子/R-S/金属盐后缀。"""
    from namepredict.layer5.stereo import apply_rs_prefix
    kind, n = _parent_n(numbered)
    names = _names_for(kind, n, numbered)
    if not names:
        return _unsupported(n, kind)

    en, zh = join_kind_name(kind, _prefix_for(numbered, kind, n), names, numbered)
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
