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
    ("anhydride", 12), ("ester", 11), ("diester", 11), ("carbamate", 11),
    ("acyl_chloride", 10), ("amide", 9),
    ("nitrile", 8),
    ("aldehyde", 7),
    ("ketone", 6), ("dione", 6), ("cycloketone", 6),
    ("alcohol", 5), ("diol", 5), ("triol", 5),
    ("cycloalcohol", 5), ("thiol", 4),
    ("amine", 3), ("diamine", 3), ("sec_amine", 3), ("tert_amine", 3),
    ("cycloamine", 3), ("phosphate", 2), ("phosphonic", 2),
    ("ether", 2), ("sulfide", 2),
)
_ARENE_NAMED: tuple[tuple[str, str, str, int], ...] = (
    ("benzoic", "benzoic acid", "苯甲酸", 13),
    ("benzoyl_chloride", "benzoyl chloride", "苯甲酰氯", 10),
    ("benzonitrile", "benzonitrile", "苯甲腈", 8),
    ("benzaldehyde", "benzaldehyde", "苯甲醛", 7),
    ("acetophenone", "acetophenone", "苯乙酮", 6),
    ("phenol", "phenol", "苯酚", 5),
    ("aniline", "aniline", "苯胺", 3),
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
_FUSED_FG: tuple[tuple[str, int, bool], ...] = (
    ("quinolinecarboxylic", 13, True), ("indolecarboxylic", 13, True),
    ("naphthalenecarboxylic", 13, True), ("indazolecarbonitrile", 8, True),
    ("indazolecarbaldehyde", 7, True), ("benzothiophenol", 5, False),
    ("quinolinol", 5, False), ("benzofuranamine", 3, False),
    ("benzothiazolamine", 3, False), ("benzimidazolamine", 3, False),
    ("benzoxazolamine", 3, False),
)
_HETERO_MONO: tuple[tuple[str, str, str], ...] = (
    ("pyridine", "pyridine", "吡啶"), ("furan", "furan", "呋喃"),
    ("thiophene", "thiophene", "噻吩"), ("pyrrole", "1H-pyrrole", "吡咯"),
    ("imidazole", "1H-imidazole", "咪唑"), ("pyrazole", "1H-pyrazole", "吡唑"),
    ("oxazole", "1,3-oxazole", "恶唑"), ("thiazole", "1,3-thiazole", "噻唑"),
    ("pyrimidine", "pyrimidine", "嘧啶"), ("pyrazine", "pyrazine", "吡嗪"),
    ("pyridazine", "pyridazine", "哒嗪"),
    ("aziridine", "aziridine", "氮杂环丙烷"),
    ("oxirane", "oxirane", "环氧乙烷"),
    ("oxolane", "oxolane", "氧杂环戊烷"),
    ("oxane", "oxane", "氧杂环己烷"),
    ("pyrrolidine", "pyrrolidine", "吡咯烷"),
    ("piperidine", "piperidine", "哌啶"),
    ("morpholine", "morpholine", "吗啉"),
    ("piperazine", "piperazine", "哌嗪"),
    ("dioxolane", "1,3-dioxolane", "1,3-二氧戊环"),
    ("dioxane", "1,4-dioxane", "1,4-二氧六环"),
    ("thiolane", "thiolane", "硫杂环戊烷"),
)
_HETERO_FUSED: tuple[tuple[str, str, str], ...] = (
    ("indole", "1H-indole", "吲哚"), ("indazole", "1H-indazole", "1H-吲唑"),
    ("benzofuran", "benzofuran", "苯并呋喃"),
    ("benzothiophene", "1-benzothiophene", "苯并[b]噻吩"),
    ("benzothiazole", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("benzoxazole", "1,3-benzoxazole", "1,3-苯并噁唑"),
    ("benzimidazole", "1H-benzimidazole", "1H-苯并咪唑"),
    ("quinoline", "quinoline", "喹啉"),
    ("isoquinoline", "isoquinoline", "异喹啉"),
    ("quinazoline", "quinazoline", "喹唑啉"),
    ("quinoxaline", "quinoxaline", "喹喔啉"),
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


def _load_fused_fg() -> None:
    for k, fg, ret in _FUSED_FG:
        _add(k, fg=fg, n=2, ret=ret)


def _load_misc_ring_fg() -> None:
    _add("cycloalkanecarboxylic", fg=13)
    _add("benzenediol", fg=5)
    _add("pyridinol", fg=5)
    _add("pyridinamine", fg=3)
    _add("pyrimidinamine", fg=3)
    _add("benzenediamine", fg=3, ring="carbo", n=1)


def _load_hetero_mono() -> None:
    for kind, en, zh in _HETERO_MONO:
        _add(kind, en=en, zh=zh, ring="hetero", n=1, ret=True)


def _load_hetero_fused() -> None:
    for kind, en, zh in _HETERO_FUSED:
        _add(kind, en=en, zh=zh, ring="hetero", n=2, ret=True)


def _load_carbo_rings() -> None:
    _add("benzene", en="benzene", zh="苯", ring="carbo", n=1, ret=True)
    _add("naphthalene", en="naphthalene", zh="萘", ring="carbo", n=2, ret=True)
    _add("anthracene", en="anthracene", zh="蒽", ring="carbo", n=3, ret=True)
    _add("cycloalkane", ring="carbo", n=1)
    _add("cycloalkene", ring="carbo", n=1)
    _add("cyclopolyene", ring="carbo", n=1)


def _load_sat_hetero_repl() -> None:
    # stems filled at runtime on parent (stem_en/stem_zh); placeholders for lint
    _add(
        "sat_hetero_repl", en="heterocycloalkane", zh="杂环烷",
        ring="hetero", n=1, ret=False,
    )


def _bootstrap() -> None:
    _load_chain_fg()
    _load_arene_fg_names()
    _load_h5_cooh()
    _load_sat_cooh()
    _load_fused_fg()
    _load_misc_ring_fg()
    _load_hetero_mono()
    _load_hetero_fused()
    _load_carbo_rings()
    _load_sat_hetero_repl()


_bootstrap()
