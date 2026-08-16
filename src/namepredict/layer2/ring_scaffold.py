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
# RDKit 从 smiles 自动算，retained=True。新增环系只在此表加一条。
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
    # saturated monohetero rings（P-22.2.2；radical/取代基须识别为环而非开链）
    "pyrrolidine": {"smiles": "C1CCNC1",  "stem_en": "pyrrolidine", "stem_zh": "吡咯烷", "naming_class": "monohetero"},
    "piperidine":  {"smiles": "C1CCNCC1", "stem_en": "piperidine",  "stem_zh": "哌啶",   "naming_class": "monohetero"},
    "morpholine":  {"smiles": "C1COCCN1", "stem_en": "morpholine",  "stem_zh": "吗啉",   "naming_class": "monohetero"},
    "piperazine":  {"smiles": "C1CNCCN1", "stem_en": "piperazine",  "stem_zh": "哌嗪",   "naming_class": "monohetero"},
    "oxolane":     {"smiles": "C1CCOC1",  "stem_en": "oxolane",     "stem_zh": "四氢呋喃", "naming_class": "monohetero"},
    "oxane":       {"smiles": "C1CCCOC1", "stem_en": "oxane",       "stem_zh": "四氢吡喃", "naming_class": "monohetero"},
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

# 保留 fused 母体的固定编号标签（P-25.4）：融合桥头用字母 locant（3a/7a、4a/8a）。
FUSED56_LABELS: tuple[str, ...] = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")
NAPH_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")

# 不对称 fused 环（含杂原子）与 1,3-二唑的标准编号：模板原子索引按固定 locant 顺序排列。
# 起点为最优先杂原子，沿环编号绕开融合桥头（桥头只得字母位）。
# 对称碳环（naphthalene/anthracene）与单杂环（pyrrole/pyridine 等单杂原子）不在此列，
# 走 P-14.4 通用枚举；1,3-二唑（咪唑/吡唑/噻唑/噁唑）IUPAC 编号固定（N 必须得 1,3 位、
# 带 H/取代基的 N 为 1），通用枚举会被 principal 最小化翻转，故登记标准编号。
# _STANDARD_ORDERS[spec] = 模板原子按 locant 顺序；_STANDARD_LABELS[spec] = 对应 locant 标签。
_STANDARD_ORDERS: dict[str, tuple[int, ...]] = {
    # naph_family（10 原子）：1,2,3,4,4a,5,6,7,8,8a
    "quinoline":    (4, 5, 6, 7, 8, 9, 0, 1, 2, 3),
    "quinazoline":  (4, 5, 6, 7, 8, 9, 0, 1, 2, 3),
    # fused56（9 原子）：1,2,3,3a,4,5,6,7,7a；杂原子(1)走远离桥头方向，苯环从 3a 起
    "indole":       (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzofuran":   (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzothiophene": (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzothiazole": (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzoxazole":  (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "indazole":     (6, 5, 4, 3, 2, 1, 0, 8, 7),
    # monohetero 1,3-二唑（5 原子）：N1(带 H/取代)/C2/N3/C4/C5。
    # 噻唑/噁唑杂原子异元素（S/O vs N）模板匹配无歧义，登记固定编号；
    # 咪唑/吡唑双 N 对称（模板 [nH] 对两个 N 可互换匹配），N1 须动态取
    # 带取代基/H 的 N，走 numbering_engine._ring_hetero_start + 杂原子最小化。
    "thiazole":     (2, 3, 4, 0, 1),
    "oxazole":      (2, 3, 4, 0, 1),
}
_STANDARD_LABELS: dict[str, tuple[str, ...]] = {
    "quinoline": NAPH_LABELS,
    "quinazoline": NAPH_LABELS,
    "indole": FUSED56_LABELS,
    "benzofuran": FUSED56_LABELS,
    "benzothiophene": FUSED56_LABELS,
    "benzothiazole": FUSED56_LABELS,
    "benzoxazole": FUSED56_LABELS,
    "indazole": FUSED56_LABELS,
    "thiazole": ("1", "2", "3", "4", "5"),
    "oxazole": ("1", "2", "3", "4", "5"),
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
    std_labels = _STANDARD_LABELS.get(sid, ())
    numbering = NumberingPolicy(
        mode=_NUMBERING_MODE.get(entry["naming_class"], "fixed"),
        standard_path=std_labels,
    )
    return ScaffoldSpec(
        id=sid, naming_class=entry["naming_class"],
        stem_en=entry["stem_en"], stem_zh=entry["stem_zh"],
        n_rings=n_rings, ring=ring, retained=True,
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
        "aromatic": all(q.GetAtomWithIdx(i).GetIsAromatic() for i in range(q.GetNumAtoms())),
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
    hit = _match_with_map(info, atom_ids)
    return hit[0] if hit else None


def _match_with_map(info: dict, atom_ids) -> tuple[str, tuple[int, ...]] | None:
    """模板精确覆盖 atom_ids 时返回 (sid, match)；match[i]=模板原子 i 对应的分子原子。

    match 供固定编号（standard_path）把模板原子映射到分子原子。
    """
    mol = info["mol"]
    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)
    for sid, q in _Q.items():
        if _TEMPLATE_ELEM[sid] != elem:
            continue
        for m in mol.GetSubstructMatches(q, uniquify=True):
            if set(m) == atoms:
                return sid, m
    return None


def standard_chain(spec_id: str | None, match: tuple[int, ...] | None) -> list[int] | None:
    """把 fused 模板固定编号映射到分子：返回分子原子按标准 locant 顺序的列表。

    无标准顺序或 match 缺失/长度不符时返回 None（走 P-14.4 通用枚举）。
    """
    if not match:
        return None
    order = _STANDARD_ORDERS.get(spec_id or "")
    if not order or len(order) != len(match):
        return None
    return [match[t] for t in order]


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
