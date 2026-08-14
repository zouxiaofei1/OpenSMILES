"""ParentKind metadata registry — scoring + L5 stem authority (P-44)."""
from __future__ import annotations

from dataclasses import dataclass, replace

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass as FG
from namepredict.layer2.principal import legacy_rank


# 数量派生 kind(diacid/diol/diamine/…)已删:multiplicity 由 principal_expression_facts 承载。
_KIND_CLASS = {
    "acid": FG.ACID,
    "ester": FG.ESTER,
    "amide": FG.AMIDE, "nitrile": FG.NITRILE, "aldehyde": FG.ALDEHYDE,
    "ketone": FG.KETONE, "dione": FG.KETONE,
    "alcohol": FG.ALCOHOL,
    "amine": FG.AMINE,
    "sec_amine": FG.AMINE, "tert_amine": FG.AMINE,
    "tetraalkylammonium": FG.QUATERNARY_AMMONIUM,
    "alkane": FG.NONE,

    #暂时无效的kind
    "ether": FG.ETHER,"thiol": FG.THIOL, "anhydride": FG.ANHYDRIDE,"acyl_chloride": FG.ACYL_HALIDE, "acyl_bromide": FG.ACYL_HALIDE,
}


def _kind_class(kind: str) -> FG | None:
    return _KIND_CLASS.get(kind)


def _principal_rank(kind: str, fallback: int = 0) -> int:
    group_class = _kind_class(kind)
    return legacy_rank(group_class) if group_class else fallback


@dataclass(frozen=True)
class KindMeta:
    kind: str
    en: str | None = None
    zh: str | None = None
    fg_rank: int = 0
    ring: str = "none"  # "none" | "hetero" | "carbo"
    n_rings: int = 0
    retained: bool = False


_REG: dict[str, KindMeta] = {}

def register(meta: KindMeta) -> None:
    rank = _principal_rank(meta.kind, meta.fg_rank)
    _REG[meta.kind] = replace(meta, fg_rank=rank)


def get(kind: str) -> KindMeta | None:
    return _REG.get(kind)


def fg_rank(kind: str) -> int:
    meta = get(kind)
    return _principal_rank(kind, meta.fg_rank if meta else 0)


def has_principal_fg(kind: str) -> int:
    return 1 if fg_rank(kind) else 0


def is_hetero_ring(kind: str) -> int:
    m = get(kind)
    return 1 if m and m.ring == "hetero" else 0


def is_carbo_ring(kind: str) -> int:
    m = get(kind)
    return 1 if m and m.ring == "carbo" else 0


def n_rings_of(kind: str) -> int:
    m = get(kind)
    return m.n_rings if m else 0


def retained_bonus(kind: str) -> int:
    m = get(kind)
    return 1 if m and m.retained else 0


def parent_names(kind: str) -> tuple[str, str] | None:
    m = get(kind)
    if m is None or m.en is None or m.zh is None:
        return None
    return m.en, m.zh


def _attach_numbering_scaffold(packed: dict) -> dict:
    from namepredict.layer2.ring_scaffold import numbering_scaffold_facts
    facts = numbering_scaffold_facts(
        packed.get("scaffold_id") or packed.get("kind"), len(packed.get("chain") or ()),
    )
    return packed if facts is None else {
        **packed, "numbering_scaffold": facts, "numbering_scaffold_required": True,
    }


def pack_parent_stem(parent: dict, mol=None) -> dict:
    packed = parent if parent.get("mol") is not None else {**parent, "mol": mol}
    names = parent_names(packed.get("kind") or "")
    if names is not None and not (packed.get("stem_en") or packed.get("stem_zh")):
        packed = {**packed, "stem_en": names[0], "stem_zh": names[1]}
    return _attach_numbering_scaffold(packed)


def all_kinds() -> frozenset[str]:
    return frozenset(_REG)


def _add(
    kind: str,
    *,
    en: str | None = None,
    zh: str | None = None,
    fg: int = 0,
    ring: str = "none",
    n: int = 0,
    ret: bool = False,
) -> None:
    register(KindMeta(kind, en, zh, fg, ring, n, ret))


def _load_chain_fg() -> None:
    for k in _KIND_CLASS:
        _add(k)


def _spec_to_meta(sp) -> KindMeta:
    return KindMeta(
        sp.id, sp.stem_en, sp.stem_zh, sp.fg_rank, sp.ring, sp.n_rings,
        sp.retained,
    )


def _load_from_scaffold_specs() -> None:
    """Single stem authority: ScaffoldSpec → KindMeta (overwrite if present)."""
    from namepredict.layer2.ring_scaffold import all_specs

    for sp in all_specs():
        if sp.stem_en is None or sp.stem_zh is None:
            continue
        register(_spec_to_meta(sp))


def _bootstrap() -> None:
    _load_chain_fg()
    _load_from_scaffold_specs()  # last: Spec is stem authority


_bootstrap()
