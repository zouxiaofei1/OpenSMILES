"""骨架规格 + 保留 SMILES 模板注册表 + 环解析（2026-08-15 三合一；_TEMPLATES 为唯一事实来源，派生全部 ScaffoldSpec/ScaffoldIdentity）。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol, MolFromSmiles

from namepredict.layer2.parent_skeleton import ParentSkeleton


# ScaffoldSpec 定义（命名/编号元数据；仅 L2 数据，不做命名组装）

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
        """返回本规格的 ScaffoldIdentity。"""
        return identity_of(self)

ScaffoldId = str


@dataclass(frozen=True)
class ScaffoldIdentity:
    id: ScaffoldId
    naming_class: str
    n_rings: int
    ring: str


def identity_of(spec) -> ScaffoldIdentity:
    """由 spec 构造 ScaffoldIdentity。"""
    return ScaffoldIdentity(spec.id, spec.naming_class, spec.n_rings, spec.ring)


# 保留母体 SMILES 模板注册表（唯一事实来源；原 specs.py + retained_templates.py 合并）

# 每个保留母体一条：smiles 模板（`match_retained` 子图同构识别）+ 命名元数据
# （词干 / 命名类）。`_spec_from_template` 派生 ScaffoldSpec，n_rings/ring 由
# RDKit 从 smiles 自动算，retained=True、fg_rank=0。新增环系只在此表加一条。
_NUMBERING_MODE = {
    "mono_carbo": "fixed_roles",
    "monohetero": "fixed_hetero",
    "fused56": "fused56_fixed",
    "naph_family": "naph_family",
    "anthra": "anthra",
}

_TEMPLATES: dict[str, dict] = {
    # carbocycles
    "benzene":     {"smiles": "c1ccccc1",             "stem_en": "benzene",    "stem_zh": "苯",   "naming_class": "mono_carbo"},
    "naphthalene": {"smiles": "c1ccc2ccccc2c1",       "stem_en": "naphthalene","stem_zh": "萘",    "naming_class": "naph_family"},
    "anthracene":  {"smiles": "c1ccc2cc3ccccc3cc2c1", "stem_en": "anthracene", "stem_zh": "蒽",    "naming_class": "anthra"},
    # monocyclic heteroarenes
    "furan":       {"smiles": "c1ccoc1",    "stem_en": "furan",       "stem_zh": "呋喃",   "naming_class": "monohetero"},
    "thiophene":   {"smiles": "c1ccsc1",    "stem_en": "thiophene",   "stem_zh": "噻吩",   "naming_class": "monohetero"},
    "pyrrole":     {"smiles": "c1cc[nH]c1", "stem_en": "pyrrole",     "stem_zh": "吡咯",   "naming_class": "monohetero"},
    "pyridine":    {"smiles": "n1ccccc1",   "stem_en": "pyridine",    "stem_zh": "吡啶",   "naming_class": "monohetero"},
    "pyridazine":  {"smiles": "c1ccnnc1",   "stem_en": "pyridazine",  "stem_zh": "哒嗪",   "naming_class": "monohetero"},
    "pyrimidine":  {"smiles": "c1cncnc1",   "stem_en": "pyrimidine",  "stem_zh": "嘧啶",   "naming_class": "monohetero"},
    "pyrazine":    {"smiles": "c1cnccn1",   "stem_en": "pyrazine",    "stem_zh": "吡嗪",   "naming_class": "monohetero"},
    "imidazole":   {"smiles": "c1cnc[nH]1", "stem_en": "imidazole",   "stem_zh": "咪唑",   "naming_class": "monohetero"},
    "pyrazole":    {"smiles": "c1ccn[nH]1", "stem_en": "pyrazole",    "stem_zh": "吡唑",   "naming_class": "monohetero"},
    "oxazole":     {"smiles": "c1cocn1",    "stem_en": "oxazole",     "stem_zh": "噁唑",   "naming_class": "monohetero"},
    "thiazole":    {"smiles": "c1cscn1",    "stem_en": "thiazole",    "stem_zh": "噻唑",   "naming_class": "monohetero"},
    # fused 5+6
    "indole":         {"smiles": "c1ccc2[nH]ccc2c1", "stem_en": "1H-indole",      "stem_zh": "吲哚",     "naming_class": "fused56"},
    "indazole":       {"smiles": "c1ccc2cn[nH]c2c1", "stem_en": "indazole",       "stem_zh": "吲唑",     "naming_class": "fused56"},
    "benzimidazole":  {"smiles": "c1ccc2[nH]cnc2c1", "stem_en": "benzimidazole",  "stem_zh": "苯并咪唑", "naming_class": "fused56"},
    "benzofuran":     {"smiles": "c1ccc2occc2c1",    "stem_en": "benzofuran",     "stem_zh": "苯并呋喃", "naming_class": "fused56"},
    "benzothiophene": {"smiles": "c1ccc2sccc2c1",    "stem_en": "benzothiophene", "stem_zh": "苯并噻吩", "naming_class": "fused56"},
    "benzothiazole":  {"smiles": "c1ccc2scnc2c1",    "stem_en": "benzothiazole",  "stem_zh": "苯并噻唑", "naming_class": "fused56"},
    "benzoxazole":    {"smiles": "c1ccc2ocnc2c1",    "stem_en": "benzoxazole",    "stem_zh": "苯并噁唑", "naming_class": "fused56"},
    # fused 6+6
    "quinoline":    {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline",    "stem_zh": "喹啉",   "naming_class": "naph_family"},
    "isoquinoline": {"smiles": "c1nccc2ccccc21", "stem_en": "isoquinoline", "stem_zh": "异喹啉", "naming_class": "naph_family"},
    "quinazoline":  {"smiles": "c1ccc2ncncc2c1", "stem_en": "quinazoline",  "stem_zh": "喹唑啉", "naming_class": "naph_family"},
    "quinoxaline":  {"smiles": "c1ccc2nccnc2c1", "stem_en": "quinoxaline",  "stem_zh": "喹噁啉", "naming_class": "naph_family"},
    # NOTE: carbonyl mothers（benzoquinone / anthraquinone / chromenone /
    # ortho_benzoquinone）不入表：模板含环外 =O，匹配集会超出环系统原子集。
}

# 查询子结构与元素签名，import 时构建一次。
_Q: dict[str, Mol] = {sid: MolFromSmiles(entry["smiles"]) for sid, entry in _TEMPLATES.items()}


def _elem_sig(mol: Mol, atom_ids) -> frozenset:
    """原子集的元素组成签名（(Z, 计数) 冻结集合）。"""
    return frozenset(Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids).items())


_TEMPLATE_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q.items()}


def _spec_from_template(sid: str, entry: dict) -> ScaffoldSpec:
    """由模板条目派生 ScaffoldSpec（n_rings/ring 从 smiles 自动算）。"""
    q = _Q[sid]
    n_rings = len(q.GetRingInfo().AtomRings())
    ring = "carbo" if all(q.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in range(q.GetNumAtoms())) else "hetero"
    numbering = NumberingPolicy(mode=_NUMBERING_MODE.get(entry["naming_class"], "fixed"))
    return ScaffoldSpec(
        id=sid, naming_class=entry["naming_class"],
        stem_en=entry["stem_en"], stem_zh=entry["stem_zh"],
        n_rings=n_rings, ring=ring, retained=True, fg_rank=0,
        numbering=numbering,
    )


_ALL_SPECS: tuple[ScaffoldSpec, ...] = tuple(
    _spec_from_template(sid, entry) for sid, entry in _TEMPLATES.items()
)
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}
_IDENTITIES: dict[str, ScaffoldIdentity] = {s.id: s.identity for s in _ALL_SPECS}


def get_identity(spec_id: str) -> ScaffoldIdentity | None:
    """按 id 查 ScaffoldIdentity（无则 None）。"""
    return _IDENTITIES.get(spec_id)


def all_identities() -> tuple[ScaffoldIdentity, ...]:
    """返回全部 ScaffoldIdentity。"""
    return tuple(_IDENTITIES.values())


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    """按 id 查 ScaffoldSpec（无则 None）。"""
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    """返回全部 ScaffoldSpec（由 _TEMPLATES 派生）。"""
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
    """使用 fused56 位次标签的 kind / scaffold id。"""
    return kind_ids_for("fused56")


def naph_kind_ids() -> frozenset[str]:
    """使用 naph 位次标签的 kind / scaffold id。"""
    return kind_ids_for("naph_family")


def monohetero_kind_ids() -> frozenset[str]:
    """单杂环保留骨架的 kind / scaffold id。"""
    return kind_ids_for("monohetero")


def kind_ids_for(naming_class: str) -> frozenset[str]:
    """给定 naming_class（fused56 / naph_family / monohetero / …）的 id 集合。"""
    return frozenset(s.id for s in _ALL_SPECS if s.naming_class == naming_class)


# RetainedEntry（保留母体条目；由 _TEMPLATES 派生，供 registry/get_entry/match_systems 查询）

RetainedEntry = dict


def _entry_for(sid: str) -> RetainedEntry | None:
    """由模板条目派生保留条目（拓扑字段 + 词干）。"""
    entry = _TEMPLATES.get(sid)
    if entry is None:
        return None
    q = _Q[sid]
    hetero_Z = tuple(sorted(
        q.GetAtomWithIdx(i).GetAtomicNum() for i in range(q.GetNumAtoms())
        if q.GetAtomWithIdx(i).GetAtomicNum() != 6
    ))
    n_rings = len(q.GetRingInfo().AtomRings())
    return {
        "kind": sid,
        "n_rings": n_rings,
        "n_atoms": q.GetNumAtoms(),
        "hetero_Z": hetero_Z,
        "topology": "fused" if n_rings > 1 else "mono",
        "aromatic": True,
        "en": entry["stem_en"],
        "zh": entry["stem_zh"],
    }


def registry() -> dict[str, RetainedEntry]:
    """返回带词干的保留拓扑注册表（由 _TEMPLATES 派生）。"""
    return {sid: _entry_for(sid) for sid in _TEMPLATES}


def get_entry(scaffold_id: str) -> RetainedEntry | None:
    """按 scaffold_id 查带词干的保留条目。"""
    return _entry_for(scaffold_id)


def match_systems(info: dict) -> list[tuple[str, dict, RetainedEntry]]:
    """返回匹配系统的 (scaffold_id, ring_system, entry)（模板子图同构语义）。"""
    out: list[tuple[str, dict, RetainedEntry]] = []
    for system in info.get("ring_systems") or []:
        sid = match_retained(info, system.get("atom_ids") or ())
        if sid:
            out.append((sid, system, _entry_for(sid)))
    return out


def match_scaffold_ids(info: dict) -> list[str]:
    """返回匹配的 scaffold id 列表。"""
    return [sid for sid, _, _ in match_systems(info)]


def match_retained(info: dict, atom_ids) -> str | None:
    """返回模板精确覆盖 atom_ids 的保留母体 sid；无命中返回 None。

    元素签名预过滤跳过组成不符的模板，再跑子图同构。多个模板同命中时
    按模板表顺序取第一个（元素标注下实际不会发生，防御性兜底）。
    """
    mol = info["mol"]
    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)
    for sid, q in _Q.items():
        if _TEMPLATE_ELEM[sid] != elem:
            continue
        if any(set(m) == atoms for m in mol.GetSubstructMatches(q, uniquify=True)):
            return sid
    return None


# 环解析

def _matched_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    """按模板子图同构匹配骨架的 scaffold id。"""
    return match_retained(info, skeleton.atom_ids)


def _generic_carbocycle(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    """纯碳环兜底为 carbocycle 身份。"""
    mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids):
        return None
    return ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")


def resolve_ring_scaffold(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    """解析骨架的 scaffold 身份（显式/模板匹配/兜底碳环）。

    模板命中即解析出该母体的 ScaffoldIdentity（_TEMPLATES 唯一来源派生）；
    无模板命中时全碳环兜底 carbocycle，杂环返回 None。
    """
    direct = get_identity(skeleton.scaffold_id or "")
    if direct:
        return direct
    sid = _matched_id(info, skeleton)
    if sid:
        identity = get_identity(sid)
        if identity:
            return identity
    return _generic_carbocycle(info, skeleton)
