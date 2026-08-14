from __future__ import annotations
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


def _benzenediol_names(numbered: dict) -> tuple[str, str] | None:
    """Benzene parent + 2 hydroxyls: systematic benzene-N,N-diol (P-22.1.3)."""
    rec = next((f for f in numbered.get("fg_locants") or [] if f.get("kind") == "oh"), None)
    locs = rec.get("locants") if rec else None
    if not locs or len(locs) != 2:
        return None
    loc_s = ",".join(str(x) for x in locs)
    return f"benzene-{loc_s}-diol", f"苯-{loc_s}-二酚"


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Chain-engine dispatch for table kinds, then specific workers."""
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
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
