from __future__ import annotations
from namepredict.layer5.chain_engine import _KIND_TABLE, _alkane_names, _chain_names
from namepredict.layer5.typed_kinds import _typed_expression_kind
from namepredict.layer5.benzene_names import join_kind_name
from namepredict.layer5.stems import maybe_anion_names, maybe_metal_salt_names
from namepredict.layer5.assembler_prefixes import _prefix_for
from namepredict.types import NameResult

def _fail(meta: dict | None = None) -> NameResult: return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})
def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult: return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)

def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Chain-engine dispatch for table kinds, then specific workers."""
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
        return _chain_names(entry, n, numbered)
  

    if kind == "phenyl":
        return "phenyl", "苯基"
    if kind == "benzene":
        return ("benzene", "苯")
    if kind == "benzoate":
        return ("benzoate", "苯甲酸")
   
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
