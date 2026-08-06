"""Scaffold specs registry (L2 data only; no naming assembly)."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2.scaffold.identity import ScaffoldIdentity, identity_of


@dataclass(frozen=True)
class NumberingPolicy:
    mode: str
    standard_path: tuple = ()
    materialize_plan: bool = True
    anchors: tuple[str, ...] = ()
    substitutable: frozenset[str] = frozenset()


@dataclass(frozen=True)
class ScaffoldSpec:
    id: str
    naming_class: str
    stem_en: str | None
    stem_zh: str | None
    n_rings: int
    ring: str
    retained: bool
    fg_rank: int
    numbering: NumberingPolicy
    sub_rules: object | None = None
    principal_slots: object | None = None

    @property
    def identity(self) -> ScaffoldIdentity:
        return identity_of(self)


# Shared fused 5+6 path labels (IUPAC P-22.2.1 / P-25): hetero=1 … 7a.
FUSED56_LABELS: tuple[str, ...] = (
    "1", "2", "3", "3a", "4", "5", "6", "7", "7a",
)
# Naphthalene / quinoline family path labels (P-25): 1…4a…8a.
NAPH_LABELS: tuple[str, ...] = (
    "1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a",
)
ANTHRA_LABELS: tuple[str, ...] = (
    "1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a",
)


def _carbo(sid: str, mode: str) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="carbocycle", stem_en=None, stem_zh=None,
        n_rings=1, ring="carbo", retained=False, fg_rank=0, numbering=pol,
    )


def _fused56(
    sid: str, stem_en: str | None, stem_zh: str | None, fg_rank: int = 0,
    *, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fused56_fixed", standard_path=FUSED56_LABELS)
    return ScaffoldSpec(
        id=sid, naming_class="fused56", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring="hetero", retained=retained, fg_rank=fg_rank,
        numbering=pol,
    )


def _naph(
    sid: str, stem_en: str | None, stem_zh: str | None, *,
    ring: str = "hetero", fg_rank: int = 0, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="naph_family", standard_path=NAPH_LABELS)
    return ScaffoldSpec(
        id=sid, naming_class="naph_family", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring=ring, retained=retained, fg_rank=fg_rank, numbering=pol,
    )


def _monohetero(
    sid: str, stem_en: str, stem_zh: str, *, fg_rank: int = 0,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fixed_hetero")
    return ScaffoldSpec(
        id=sid, naming_class="monohetero", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="hetero", retained=True, fg_rank=fg_rank,
        numbering=pol,
    )


def _mono_carbo(
    sid: str, stem_en: str, stem_zh: str, *, mode: str = "fixed_roles",
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="mono_carbo", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="carbo", retained=True, fg_rank=0, numbering=pol,
    )


def _poly_carbo(
    sid: str, stem_en: str, stem_zh: str, n_rings: int, mode: str,
    *, fg_rank: int = 0,
) -> ScaffoldSpec:
    pol = NumberingPolicy(
        mode=mode, standard_path=ANTHRA_LABELS if mode == "anthracene_fixed" else (),
    )
    return ScaffoldSpec(
        id=sid, naming_class="poly_carbo", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=n_rings, ring="carbo", retained=True, fg_rank=fg_rank,
        numbering=pol,
    )


def _benzodiazine(
    sid: str, stem_en: str, stem_zh: str, *, fg_rank: int = 0,
) -> ScaffoldSpec:
    """6+6 diazine labels are specified but no producer emits a plan yet."""
    pol = NumberingPolicy(
        mode="naph_family", standard_path=NAPH_LABELS, materialize_plan=False,
    )
    return ScaffoldSpec(
        id=sid, naming_class="benzodiazine", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring="hetero", retained=True, fg_rank=fg_rank, numbering=pol,
    )


CARBOCYCLE_SPECS: tuple[ScaffoldSpec, ...] = (
    _carbo("cycloalkane", "carbocycle_free"),
    _carbo("cycloalkene", "carbocycle_free"),
    _carbo("cyclopolyene", "poly_unsat"),
    ScaffoldSpec("cycloalkane_polycarboxylic", "carbocycle", "cycloalkane", "环烷烃", 1, "carbo", False, 13, NumberingPolicy("carbocycle_free")),
)

# O/S 5+6 retained (prior Phase 2.1) + aza 5+6 (indole / indazole / bim).
FUSED56_SPECS: tuple[ScaffoldSpec, ...] = (
    _fused56("benzofuran", "benzofuran", "苯并呋喃"),
    _fused56("benzofuranamine", "benzofuranamine", "苯并呋喃胺", fg_rank=3),
    _fused56("benzothiophene", "1-benzothiophene", "苯并[b]噻吩"),
    _fused56("benzothiophenol", "1-benzothiophenol", "苯并[b]噻吩酚", fg_rank=5),
    _fused56("benzothiazole", "1,3-benzothiazole", "1,3-苯并噻唑"),
    _fused56("benzothiazolamine", "benzothiazolamine", "苯并噻唑胺", fg_rank=3),
    _fused56("benzoxazole", "1,3-benzoxazole", "1,3-苯并噁唑"),
    _fused56("benzoxazolamine", "benzoxazolamine", "苯并噁唑胺", fg_rank=3),
    _fused56("indole", "1H-indole", "吲哚"),
    _fused56("indolecarboxylic", "indolecarboxylic", "吲哚羧酸", fg_rank=13),
    _fused56("indazole", "1H-indazole", "1H-吲唑"),
    _fused56("indazolecarbonitrile", "indazolecarbonitrile", "吲唑甲腈", fg_rank=8),
    _fused56("indazolecarbaldehyde", "indazolecarbaldehyde", "吲唑甲醛", fg_rank=7),
    _fused56("benzimidazole", "1H-benzimidazole", "1H-苯并咪唑"),
    _fused56(
        "benzimidazolamine", "benzimidazolamine", "苯并咪唑胺", fg_rank=3,
    ),
)

# Quinoline / isoquinoline / naphthalene / chromen-2-one (10-atom path).
NAPH_FAMILY_SPECS: tuple[ScaffoldSpec, ...] = (
    _naph("quinoline", "quinoline", "喹啉"),
    _naph("isoquinoline", "isoquinoline", "异喹啉"),
    _naph("quinolinol", "quinolinol", "喹啉酚", fg_rank=5),
    _naph(
        "quinolinecarboxylic", "quinolinecarboxylic", "喹啉羧酸", fg_rank=13,
    ),
    _naph("naphthalene", "naphthalene", "萘", ring="carbo"),
    _naph(
        "naphthalenecarboxylic", "naphthalenecarboxylic", "萘羧酸",
        ring="carbo", fg_rank=13,
    ),
    # Arene FG parents: OH/NH2 on fused carbo/hetero rings (P-63.1.4 / P-62.2.1).
    _naph("naphthalenol", "naphthalenol", "萘酚", ring="carbo", fg_rank=5),
    _naph("naphthalenediol", "naphthalenediol", "萘二酚", ring="carbo", fg_rank=5),
    _naph("naphthalenamine", "naphthalenamine", "萘胺", ring="carbo", fg_rank=3),
    _naph("quinolinediol", "quinolinediol", "喹啉二酚", fg_rank=5),
    # Arene FG parents: CHO/CN on fused rings (P-66.6.1 / P-66.5.1).
    _naph("naphthalenecarbaldehyde", "naphthalenecarbaldehyde", "萘甲醛", ring="carbo", fg_rank=7),
    _naph("quinolinecarbaldehyde", "quinolinecarbaldehyde", "喹啉甲醛", fg_rank=7),
    _naph("naphthalenecarbonitrile", "naphthalenecarbonitrile", "萘甲腈", ring="carbo", fg_rank=8),
    _naph("quinolinecarbonitrile", "quinolinecarbonitrile", "喹啉甲腈", fg_rank=8),
    # Coumarin lactone retained: EN chromen-2-one; ZH 香豆素 (ketone-class rank).
    _naph("chromenone", "chromen-2-one", "香豆素", fg_rank=6),
)

# 6+6 benzodiazines: same 10-atom naph labels; not yet in L4 Q_KINDS.
BENZODIAZINE_SPECS: tuple[ScaffoldSpec, ...] = (
    _benzodiazine("quinazoline", "quinazoline", "喹唑啉"),
    _benzodiazine("quinazolinamine", "quinazolinamine", "喹唑啉胺", fg_rank=3),
    _benzodiazine("quinoxaline", "quinoxaline", "喹喔啉"),
)

# Mono-hetero retained + Hantzsch–Widman (stem authority for kind_registry).
MONO_HETERO_SPECS: tuple[ScaffoldSpec, ...] = (
    _monohetero("pyridine", "pyridine", "吡啶"),
    _monohetero("furan", "furan", "呋喃"),
    _monohetero("thiophene", "thiophene", "噻吩"),
    _monohetero("pyrrole", "1H-pyrrole", "吡咯"),
    _monohetero("imidazole", "1H-imidazole", "咪唑"),
    _monohetero("pyrazole", "1H-pyrazole", "吡唑"),
    _monohetero("pyrazolamine", "pyrazolamine", "吡唑胺", fg_rank=3),
    _monohetero("oxazole", "1,3-oxazole", "恶唑"),
    _monohetero("thiazole", "1,3-thiazole", "噻唑"),
    _monohetero("thiazolamine", "thiazolamine", "噻唑胺", fg_rank=3),
    _monohetero("pyrimidine", "pyrimidine", "嘧啶"),
    _monohetero("pyrazine", "pyrazine", "吡嗪"),
    _monohetero("pyridazine", "pyridazine", "哒嗪"),
    _monohetero("aziridine", "aziridine", "氮杂环丙烷"),
    _monohetero("oxirane", "oxirane", "环氧乙烷"),
    _monohetero("oxolane", "oxolane", "氧杂环戊烷"),
    _monohetero("oxane", "oxane", "氧杂环己烷"),
    _monohetero("pyrrolidine", "pyrrolidine", "吡咯烷"),
    _monohetero("piperidine", "piperidine", "哌啶"),
    _monohetero("morpholine", "morpholine", "吗啉"),
    _monohetero("piperazine", "piperazine", "哌嗪"),
    _monohetero("dioxolane", "1,3-dioxolane", "1,3-二氧戊环"),
    _monohetero("dioxane", "1,4-dioxane", "1,4-二氧六环"),
    _monohetero("thiolane", "thiolane", "硫杂环戊烷"),
)

def _mono_carbo_fg(
    sid: str, stem_en: str, stem_zh: str, *, fg_rank: int, retained: bool = False,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fixed_roles")
    return ScaffoldSpec(
        id=sid, naming_class="mono_carbo", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="carbo", retained=retained, fg_rank=fg_rank,
        numbering=pol,
    )


# Mono / linear polycyclic carbo retained (benzene, anthracene).
MONO_CARBO_SPECS: tuple[ScaffoldSpec, ...] = (
    _mono_carbo("benzene", "benzene", "苯"),
    _mono_carbo_fg(
        "benzoquinone", "cyclohexa-2,5-diene-1,4-dione", "环己-2,5-二烯-1,4-二酮",
        fg_rank=6, retained=False,
    ),
    _mono_carbo_fg(
        "ortho_benzoquinone", "cyclohexa-3,5-diene-1,2-dione", "环己-3,5-二烯-1,2-二酮",
        fg_rank=6, retained=False,
    ),
)

POLY_CARBO_SPECS: tuple[ScaffoldSpec, ...] = (
    _poly_carbo("anthracene", "anthracene", "蒽", 3, "anthracene_fixed"),
    _poly_carbo(
        "anthraquinone", "9,10-anthraquinone", "蒽醌", 3, "anthracene_fixed",
        fg_rank=6,
    ),
)

_ALL_SPECS: tuple[ScaffoldSpec, ...] = (
    CARBOCYCLE_SPECS
    + FUSED56_SPECS
    + NAPH_FAMILY_SPECS
    + BENZODIAZINE_SPECS
    + MONO_HETERO_SPECS
    + MONO_CARBO_SPECS
    + POLY_CARBO_SPECS
)
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}
_IDENTITIES: dict[str, ScaffoldIdentity] = {s.id: s.identity for s in _ALL_SPECS}


_CYCLOALKANE_STEMS = {
    3: ("cyclopropane", "环丙烷"), 4: ("cyclobutane", "环丁烷"),
    5: ("cyclopentane", "环戊烷"), 6: ("cyclohexane", "环己烷"),
    7: ("cycloheptane", "环庚烷"), 8: ("cyclooctane", "环辛烷"),
    9: ("cyclononane", "环壬烷"), 10: ("cyclodecane", "环癸烷"),
}


def cycloalkane_polyacid_stem(n: int) -> tuple[str, str] | None:
    """Parameterized C3–C10 stem facts for the cyclo polyacid ScaffoldSpec."""
    return _CYCLOALKANE_STEMS.get(n)


def get_identity(spec_id: str) -> ScaffoldIdentity | None:
    return _IDENTITIES.get(spec_id)


def all_identities() -> tuple[ScaffoldIdentity, ...]:
    return tuple(_IDENTITIES.values())


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    return _BY_ID.get(spec_id)



def all_specs() -> tuple[ScaffoldSpec, ...]:
    return _ALL_SPECS


def numbering_scaffold_facts(spec_id: str | None, atom_count: int) -> dict | None:
    """Materialize pure parent facts; ``relative_stereo`` reserves ring-face constraints."""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.numbering.materialize_plan:
        return None
    labels = spec.numbering.standard_path
    if spec.id in {"cycloalkane", "cycloalkane_polycarboxylic"}:
        labels = tuple(str(i) for i in range(1, atom_count + 1)) if 3 <= atom_count <= 10 else ()
    return None if not labels or len(labels) != atom_count else {
        "scaffold_id": spec.id, "labels": labels, "relative_stereo": None,
    }


def fused56_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids that use fused56 1…7a labels."""
    return frozenset(s.id for s in FUSED56_SPECS)


def naph_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids that use naph 1…8a labels."""
    return frozenset(s.id for s in NAPH_FAMILY_SPECS)


def monohetero_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids for mono-hetero retained scaffolds."""
    return frozenset(s.id for s in MONO_HETERO_SPECS)


def kind_ids_for(naming_class: str) -> frozenset[str]:
    """Ids for a naming_class (fused56 / naph_family / monohetero / …)."""
    return frozenset(s.id for s in _ALL_SPECS if s.naming_class == naming_class)
