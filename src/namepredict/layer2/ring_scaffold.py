"""骨架规格 + 保留拓扑注册表 + 环解析（2026-08-14 由 ``specs.py`` / ``retained_registry.py`` / ``ring_scaffold.py`` 合并；ScaffoldSpec/NumberingPolicy 是保留骨架词干与编号策略的单一来源）。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2.identity import ScaffoldIdentity, identity_of
from namepredict.layer2.parent_skeleton import ParentSkeleton


# ScaffoldSpec 注册表（原 specs.py；仅 L2 数据，不做命名组装）

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


# 共用的稠合 5+6 位次标签（IUPAC P-22.2.1 / P-25）：hetero=1 … 7a。
FUSED56_LABELS: tuple[str, ...] = (
    # 位次标签："1", "2", "3", "3a", "4", "5", "6", "7", "7a",
)
# 萘 / 喹啉家族位次标签（P-25）：1…4a…8a。
NAPH_LABELS: tuple[str, ...] = (
    # 位次标签："1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a",
)
ANTHRA_LABELS: tuple[str, ...] = (
    # 位次标签："1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a",
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
    sid: str, stem_en: str, stem_zh: str, *,
    mode: str = "fixed_roles", fg_rank: int = 0, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="mono_carbo", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="carbo", retained=retained, fg_rank=fg_rank, numbering=pol,
    )




_ALL_SPECS: tuple[ScaffoldSpec, ...] = (
    _mono_carbo("benzene", "benzene", "苯"),
    _monohetero("pyridine", "pyridine", "吡啶"),
    _naph("naphthalene", "naphthalene", "萘", ring="carbo"),
    _fused56("indole", "1H-indole", "吲哚"),
)
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}
_IDENTITIES: dict[str, ScaffoldIdentity] = {s.id: s.identity for s in _ALL_SPECS}


def get_identity(spec_id: str) -> ScaffoldIdentity | None:
    return _IDENTITIES.get(spec_id)


def all_identities() -> tuple[ScaffoldIdentity, ...]:
    return tuple(_IDENTITIES.values())


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    return _ALL_SPECS


def numbering_scaffold_facts(spec_id: str | None, atom_count: int) -> dict | None:
    """物化纯母体 facts；``relative_stereo`` 预留环面约束。"""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.numbering.materialize_plan:
        return None
    labels = spec.numbering.standard_path
    return None if not labels or len(labels) != atom_count else {
        "scaffold_id": spec.id, "labels": labels, "relative_stereo": None,
    }


def fused56_kind_ids() -> frozenset[str]:
    """使用 fused56 1…7a 位次标签的 kind / scaffold id。"""
    return kind_ids_for("fused56")


def naph_kind_ids() -> frozenset[str]:
    """使用 naph 1…8a 位次标签的 kind / scaffold id。"""
    return kind_ids_for("naph_family")


def monohetero_kind_ids() -> frozenset[str]:
    """单杂环保留骨架的 kind / scaffold id。"""
    return kind_ids_for("monohetero")


def kind_ids_for(naming_class: str) -> frozenset[str]:
    """给定 naming_class（fused56 / naph_family / monohetero / …）的 id 集合。"""
    return frozenset(s.id for s in _ALL_SPECS if s.naming_class == naming_class)


# Retained topology 注册表（原 retained_registry.py；P-22 / P-25）

# scaffold_id → topology 条目（en/zh 读取时经 get_spec 解析）；键：kind, n_rings, n_atoms, hetero_Z(排序), topology, aromatic
RetainedEntry = dict

# 仅含拓扑的表；id 必须 ⊆ ScaffoldSpec 注册表。
_TOPOLOGY: dict[str, dict] = {
    "benzene": {
        "kind": "benzene",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (),
        "topology": "mono",
        "aromatic": True,
    },
    "pyridine": {
        "kind": "pyridine",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (7,),
        "topology": "mono",
        "aromatic": True,
    },
    "naphthalene": {
        "kind": "naphthalene",
        "n_rings": 2,
        "n_atoms": 10,
        "hetero_Z": (),
        "topology": "fused",
        "aromatic": True,
    },
    "indole": {
        "kind": "indole",
        "n_rings": 2,
        "n_atoms": 9,
        "hetero_Z": (7,),
        "topology": "fused",
        "aromatic": True,
    },
}


def _entry_with_stems(sid: str, topo: dict) -> RetainedEntry:
    sp = get_spec(sid)
    en = sp.stem_en if sp else None
    zh = sp.stem_zh if sp else None
    return {**topo, "en": en, "zh": zh}


def registry() -> dict[str, RetainedEntry]:
    return {sid: _entry_with_stems(sid, t) for sid, t in _TOPOLOGY.items()}


def get_entry(scaffold_id: str) -> RetainedEntry | None:
    topo = _TOPOLOGY.get(scaffold_id)
    return None if topo is None else _entry_with_stems(scaffold_id, topo)


def _hetero_Z_tuple(system: dict) -> tuple[int, ...]:
    return tuple(sorted(h["Z"] for h in system.get("hetero_atoms") or []))


def _system_matches(system: dict, entry: RetainedEntry) -> bool:
    if system.get("n_rings") != entry["n_rings"]:
        return False
    if system.get("n_atoms") != entry["n_atoms"]:
        return False
    if system.get("topology") != entry["topology"]:
        return False
    if bool(system.get("is_aromatic_mancude")) != bool(entry.get("aromatic")):
        return False
    return _hetero_Z_tuple(system) == tuple(entry["hetero_Z"])


def match_systems(info: dict) -> list[tuple[str, dict, RetainedEntry]]:
    """返回匹配系统的 (scaffold_id, ring_system, entry)。"""
    out: list[tuple[str, dict, RetainedEntry]] = []
    for system in info.get("ring_systems") or []:
        for sid, topo in _TOPOLOGY.items():
            entry = _entry_with_stems(sid, topo)
            if _system_matches(system, entry):
                out.append((sid, system, entry))
    return out


def match_scaffold_ids(info: dict) -> list[str]:
    return [sid for sid, _, _ in match_systems(info)]


# 环解析（原 ring_scaffold.py）

def _matched_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    atoms = set(skeleton.atom_ids)
    return next((sid for sid, system, _ in match_systems(info)
                 if set(system.get("atom_ids") or ()) == atoms), None)


def _generic_carbocycle(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids):
        return None
    return ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")


def resolve_ring_scaffold(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    direct = get_identity(skeleton.scaffold_id or "")
    if direct:
        return direct
    sid = _matched_id(info, skeleton)
    if sid:
        return get_identity(sid)
    return _generic_carbocycle(info, skeleton)
