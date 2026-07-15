"""ParentKind metadata registry — scoring + L5 stem authority (P-44)."""
from __future__ import annotations

from dataclasses import dataclass


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

# --- bootstrap tables (module-level; kept out of function bodies) ---
_CHAIN_FG: tuple[tuple[str, int], ...] = (
    ("acid", 13), ("diacid", 13),
    ("anhydride", 12), ("ester", 11), ("diester", 11), ("carbamate", 11), ("carbonate", 11),
    ("sulfonate", 11),
    ("acyl_chloride", 10), ("acyl_bromide", 10), ("sulfonyl_chloride", 10),
    ("amide", 9), ("urea", 9), ("guanidine", 9), ("sulfonamide", 9), ("nitrile", 8),
    ("aldehyde", 7),
    ("ketone", 6), ("dione", 6), ("cycloketone", 6),
    ("alcohol", 5), ("diol", 5), ("triol", 5),
    ("cycloalcohol", 5), ("thiol", 4),
    ("amine", 3), ("diamine", 3), ("sec_amine", 3), ("tert_amine", 3),
    ("cycloamine", 3), ("hydrazine", 4), ("phosphate", 2), ("phosphonic", 2),
    ("ether", 2), ("sulfide", 2), ("sulfoxide", 2),
    ("isocyanate", 8), ("isothiocyanate", 8),
)
_ARENE_NAMED: tuple[tuple[str, str, str, int], ...] = (
    ("benzoic", "benzoic acid", "苯甲酸", 13),
    ("benzoyl_chloride", "benzoyl chloride", "苯甲酰氯", 10),
    ("benzoyl_bromide", "benzoyl bromide", "苯甲酰溴", 10),
    ("benzonitrile", "benzonitrile", "苯甲腈", 8),
    ("benzaldehyde", "benzaldehyde", "苯甲醛", 7),
    ("acetophenone", "acetophenone", "苯乙酮", 6),
    ("phenol", "phenol", "苯酚", 5),
    ("aniline", "aniline", "苯胺", 3),
    ("boronic", "phenylboronic acid", "苯基硼酸", 12),
)
_H5_COOH = (
    "furancarboxylic", "thiophenecarboxylic", "pyrrolecarboxylic",
    "imidazolecarboxylic", "pyrazolecarboxylic", "pyridinecarboxylic",
    "pyridinecarbonitrile",
)
_SAT_COOH = (
    "piperidinecarboxylic", "pyrrolidinecarboxylic",
    "piperazinecarboxylic", "morpholinecarboxylic",
    "oxolanecarboxylic", "oxanecarboxylic",
    "thiolanecarboxylic", "aziridinecarboxylic",
)
# Fused FG kinds without Spec stems still need fg_rank registration.
# Spec is authority when present; this table only covers residual FG tags.
_MISC_RING_FG: tuple[tuple[str, int, str, int, bool], ...] = (
    ("cycloalkanecarboxylic", 13, "none", 0, False),
    ("benzenediol", 5, "none", 0, False),
    ("pyridinol", 5, "none", 0, False),
    ("pyridinamine", 3, "none", 0, False),
    ("pyrimidinamine", 3, "none", 0, False),
    ("benzenediamine", 3, "carbo", 1, False),
)


def register(meta: KindMeta) -> None:
    _REG[meta.kind] = meta


def get(kind: str) -> KindMeta | None:
    return _REG.get(kind)


def fg_rank(kind: str) -> int:
    m = get(kind)
    return m.fg_rank if m else 0


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


def all_kinds() -> frozenset[str]:
    return frozenset(_REG)


# Ordered ring parent producers (filled by layer2.ring_producers bootstrap).
_RING_TRY: list = []
_RING_BOOTSTRAPPED = False


def register_ring_try(fn) -> None:
    """Append a ring parent producer: (info) -> parent dict | None."""
    _RING_TRY.append(fn)


def _ensure_ring_producers() -> None:
    """Import ring_producers once so ring_try_fns works without candidates import."""
    global _RING_BOOTSTRAPPED
    if _RING_BOOTSTRAPPED:
        return
    from namepredict.layer2 import ring_producers as _rp  # noqa: F401
    _RING_BOOTSTRAPPED = True


def ring_try_fns() -> list:
    """Return registered ring parent try callables (order = try order)."""
    _ensure_ring_producers()
    return list(_RING_TRY)


# Ordered FG parent producers (filled by layer2.fg_producers bootstrap).
_FG_TRY: list = []
_FG_BOOTSTRAPPED = False


def register_fg_try(fn) -> None:
    """Append an FG parent producer: (info) -> parent dict | None."""
    _FG_TRY.append(fn)


def _ensure_fg_producers() -> None:
    """Import fg_producers once so fg_try_fns works without candidates import."""
    global _FG_BOOTSTRAPPED
    if _FG_BOOTSTRAPPED:
        return
    from namepredict.layer2 import fg_producers as _fp  # noqa: F401
    _FG_BOOTSTRAPPED = True


def fg_try_fns() -> list:
    """Return registered FG parent try callables (order = try order)."""
    _ensure_fg_producers()
    return list(_FG_TRY)


# Ordered unsat parent producers (filled by layer2.unsat_producers bootstrap).
_UNSAT_TRY: list = []
_UNSAT_BOOTSTRAPPED = False


def register_unsat_try(fn) -> None:
    """Append an unsat parent producer: (info) -> parent dict | None."""
    _UNSAT_TRY.append(fn)


def _ensure_unsat_producers() -> None:
    """Import unsat_producers once so unsat_try_fns works without candidates import."""
    global _UNSAT_BOOTSTRAPPED
    if _UNSAT_BOOTSTRAPPED:
        return
    from namepredict.layer2 import unsat_producers as _up  # noqa: F401
    _UNSAT_BOOTSTRAPPED = True


def unsat_try_fns() -> list:
    """Return registered unsat parent try callables (order = try order)."""
    _ensure_unsat_producers()
    return list(_UNSAT_TRY)


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
    for k, r in _CHAIN_FG:
        _add(k, fg=r)


def _load_arene_fg_names() -> None:
    for kind, en, zh, fg in _ARENE_NAMED:
        _add(kind, en=en, zh=zh, fg=fg, ret=True)
    _add("benzoate", fg=11, ret=True)


def _load_h5_cooh() -> None:
    for k in _H5_COOH:
        fg = 8 if k.endswith("carbonitrile") else 13
        _add(k, fg=fg, ret=True)


def _load_sat_cooh() -> None:
    for k in _SAT_COOH:
        _add(k, fg=13, ret=True)


def _load_misc_ring_fg() -> None:
    for k, fg, ring, n, ret in _MISC_RING_FG:
        _add(k, fg=fg, ring=ring, n=n, ret=ret)


def _load_cyclo_rings() -> None:
    """Cyclo* kinds without stems (Spec stem is None; still need ring meta)."""
    _add("cycloalkane", ring="carbo", n=1)
    _add("cycloalkene", ring="carbo", n=1)
    _add("cyclopolyene", ring="carbo", n=1)


def _load_sat_hetero_repl() -> None:
    # stems filled at runtime on parent (stem_en/stem_zh); placeholders for lint
    _add(
        "sat_hetero_repl", en="heterocycloalkane", zh="杂环烷",
        ring="hetero", n=1, ret=False,
    )


def _spec_to_meta(sp) -> KindMeta:
    return KindMeta(
        sp.id, sp.stem_en, sp.stem_zh, sp.fg_rank, sp.ring, sp.n_rings,
        sp.retained,
    )


def _load_from_scaffold_specs() -> None:
    """Single stem authority: ScaffoldSpec → KindMeta (overwrite if present)."""
    from namepredict.layer2.scaffold.specs import all_specs

    for sp in all_specs():
        if sp.stem_en is None or sp.stem_zh is None:
            continue
        register(_spec_to_meta(sp))


def _bootstrap() -> None:
    _load_chain_fg()
    _load_arene_fg_names()
    _load_h5_cooh()
    _load_sat_cooh()
    _load_misc_ring_fg()
    _load_cyclo_rings()
    _load_sat_hetero_repl()
    _load_from_scaffold_specs()  # last: Spec is stem authority


_bootstrap()
