from __future__ import annotations
from dataclasses import replace

from namepredict.layer5.chain_engine import _KIND_TABLE, _alkane_names, _chain_names
from namepredict.layer5.typed_kinds import _typed_expression_kind
from namepredict.layer5.benzene_names import join_kind_name
from namepredict.layer5.stems import maybe_anion_names, maybe_metal_salt_names
from namepredict.layer5.assembler_prefixes import _prefix_for
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult: return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})
def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult: return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

# 苯系保留名 stem（原 L2 _ARENE_NAMED 迁此；L2 只产生 kind='benzene'）。
_ARENE_RETAINED_NAMES = {
    "benzoic": ("benzoic acid", "苯甲酸"),
    "benzaldehyde": ("benzaldehyde", "苯甲醛"),
    "benzonitrile": ("benzonitrile", "苯甲腈"),
    "benzamide": ("benzamide", "苯甲酰胺"),
    "benzoate": ("benzoate", "苯甲酸"),
    "aniline": ("aniline", "苯胺"),
    "phenol": ("phenol", "苯酚"),
}

# 稠环/杂环 scaffold 词干：+FG 后缀时覆盖链状烷烃词干（正交化，不枚举萘醇/吲哚醇组合）。
_RING_STEM = {
    "naphthalene": ("naphthalen", "萘"),
    "indole": ("indol", "吲哚"),
    "pyridine": ("pyridin", "吡啶"),
    "quinoline": ("quinolin", "喹啉"),
}

# 稠环/杂环完整 base 名（-carboxylic acid 用完整词干，如 naphthalene-1-carboxylic acid）。
_RING_BASE = {
    "naphthalene": ("naphthalene", "萘"),
    "quinoline": ("quinoline", "喹啉"),
    "pyridine": ("pyridine", "吡啶"),
}


def _benzenediol_names(numbered: dict) -> tuple[str, str] | None:
    """Benzene parent + 2 hydroxyls: systematic benzene-N,N-diol (P-22.1.3)."""
    rec = next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == "oh"), None)
    locs = rec.get("locants") if rec else None
    if not locs or len(locs) != 2:
        return None
    loc_s = ",".join(str(x) for x in locs)
    return f"benzene-{loc_s}-diol", f"苯-{loc_s}-二酚"


def _scaffold_id(numbered: dict) -> str | None:
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
    """Chain-engine dispatch for table kinds, then specific workers."""
    if kind == "acid":
        exo = _exocyclic_acid_names(n, numbered)
        if exo:
            return exo
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
        # 命名 kind 正交化：环醇/酮/胺/酸与链状共用同一 _Chain，scaffold 信息运行时注入：
        #   carbocycle → cyclo 前缀；稠环/杂环（naphthalene/indole/...）→ scaffold 词干。
        # 不再预枚举 cycloalcohol/naphthalenol/indolol 组合 kind。
        sid = _scaffold_id(numbered)
        if sid in _RING_STEM:
            # 环式 FG 的 locant omit 由 L4 算出的 omit 标志决定。
            entry = replace(entry, stem=_RING_STEM[sid], omit_rule=lambda n, loc, omit: bool(omit))
        elif sid == "carbocycle":
            entry = replace(entry, cyclic=True, omit_rule=lambda n, loc, omit: bool(omit))
        return _chain_names(entry, n, numbered)


    if kind == "phenyl":
        return "phenyl", "苯基"
    if kind == "benzene":
        return ("benzene", "苯")
    if kind in _ARENE_RETAINED_NAMES:
        return _ARENE_RETAINED_NAMES[kind]
    if kind == "benzenediol":
        return _benzenediol_names(numbered)

    stem = _parent_stem_names(numbered)
    return stem

def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None

def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)
def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})


def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    from namepredict.layer5.stereo import apply_rs_prefix
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
