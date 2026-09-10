"""为选中的骨架标注主基团表达式 facts。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum

from namepredict.layer1 import fg_registry as _fg_reg
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, _anchors
from namepredict.layer2.principal import PrincipalGroupSelection, feature_spec



def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    """构造含 chain/n_carbons/kind 的母体 dict。"""
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, **kw}


class PrincipalRelation(str, Enum):
    """主基团相对骨架的位置关系（骨架内/环外）。"""
    IN_SKELETON = "in_skeleton"
    EXOCYCLIC = "exocyclic"


class PrincipalChargeState(str, Enum):
    """主基团的电荷状态（中性/全阴离子/阴离子混合）。"""
    NEUTRAL = "neutral"
    ANION = "anion"
    MIXED = "mixed"


@dataclass(frozen=True)
class PrincipalExpressionFacts:
    """主基团表达事实：类别、个数、与骨架关系、特征/附着原子集与电荷态。"""
    group_class: FunctionalGroupClass
    multiplicity: int
    relation: PrincipalRelation
    occurrence_ids: tuple[str, ...]
    characteristic_atoms: frozenset[int]
    attachment_atoms: frozenset[int]
    charge_state: PrincipalChargeState


# 链式主官能团：kind 恒为 FG 类别名；acid/alcohol/amine/ketone 任意数量恒用基团名，其余链 FG 仅单基。
# 链/数量集合由 fg_registry 的 chain/multi 标志派生（唯一事实来源）。
_CHAIN_FG = frozenset(FunctionalGroupClass(v) for v in _fg_reg.chain_fgs())
_MULTI_FG = frozenset(FunctionalGroupClass(v) for v in _fg_reg.multi_fgs())
def _anchor_fields(group_class: FunctionalGroupClass) -> tuple[str, str] | None:
    """取基团类的 (单, 复数) anchor 字段名。"""
    if group_class is FunctionalGroupClass.NONE:
        return "none_c_idx", "none_c_idxs"
    spec = feature_spec(group_class)
    return spec.anchor_fields if spec else None


def _covered(selection: PrincipalGroupSelection, skeleton: ParentSkeleton):
    """取骨架覆盖的 occurrence 子集。"""
    ids = skeleton.covered_principal_ids
    return tuple(o for o in selection.occurrences if o.id in ids)


def _chain_kind(group_class: FunctionalGroupClass, count: int) -> str | None:
    """按链 FG 类别与个数决定 kind（不支持时 None）。"""
    if group_class is FunctionalGroupClass.NONE:
        return group_class.value if count == 0 else None
    if group_class is FunctionalGroupClass.ACYL:
        return "acyl"  # 酰基残基：羰基头为 locant 1，L5 拼 -oyl/酰（P-65.1.7.2）
    if group_class is FunctionalGroupClass.RADICAL:
        return "radical"  # 自由基连接点位次由 L4 radical_c_idx 承载，kind 恒 "radical"（L5 worker 拼 -yl）
    if group_class not in _CHAIN_FG:
        return None
    if group_class not in _MULTI_FG and count != 1:
        return None
    return group_class.value if count >= 1 else None


def _principal_fields(group_class: FunctionalGroupClass, anchors: list[int]) -> dict:
    """按锚点数填单/复数 anchor 字段。"""
    single, plural = _anchor_fields(group_class)
    return {single: anchors[0], plural: anchors} if len(anchors) == 1 else {plural: anchors}


def _expression_flags(selection: PrincipalGroupSelection, occurrences) -> dict:
    """主基团表达标志（酸全阴离子时加 anion）。"""
    if selection.group_class is not FunctionalGroupClass.ACID:
        return {}
    return {"anion": True} if occurrences and all(o.payload.get("anion") for o in occurrences) else {}


def _charge_state(occurrences) -> PrincipalChargeState:
    """按 occurrence 阴离子情况推断电荷状态。"""
    charges = [bool(o.payload.get("anion")) for o in occurrences]
    if charges and all(charges):
        return PrincipalChargeState.ANION
    return PrincipalChargeState.MIXED if any(charges) else PrincipalChargeState.NEUTRAL


def _skeletal_attachments(mol, skeleton, occurrences) -> frozenset[int]:
    """计算主官能团在骨架内的附着原子（骨架外则取邻居）。"""
    atoms = set(skeleton.atom_ids)
    anchors = {i for o in occurrences for i in o.parent_anchors}
    included = anchors & atoms
    if included or mol is None:
        return frozenset(included or anchors)
    return frozenset(n.GetIdx() for i in anchors for n in mol.GetAtomWithIdx(i).GetNeighbors()
                     if n.GetIdx() in atoms)


def _facts(selection, skeleton, occurrences, mol=None) -> PrincipalExpressionFacts:
    """汇总主基团表达式 facts（关系/附着/电荷）。"""
    characteristic = frozenset(i for o in occurrences for i in o.characteristic_atoms)
    relation = PrincipalRelation.IN_SKELETON if characteristic & set(skeleton.atom_ids) else PrincipalRelation.EXOCYCLIC
    attachment = _skeletal_attachments(mol, skeleton, occurrences)
    return PrincipalExpressionFacts(selection.group_class, len(occurrences), relation,
                                    tuple(o.id for o in occurrences), characteristic, attachment,
                                    _charge_state(occurrences))


def _parent_dict(kind: str, skeleton: ParentSkeleton, occurrences, fields: dict,
                 facts: PrincipalExpressionFacts) -> dict:
    """按骨架构造带表达式 facts 的母体 dict。"""
    return {"kind": kind, "chain": list(skeleton.atom_ids), "n_carbons": len(skeleton.atom_ids),
            "covered_principal_ids": tuple(o.id for o in occurrences),
            "principal_group_count": len(occurrences), "principal_expression_facts": facts, **fields}


# 苯系保留名已迁往 L5 chain_engine 苯 variant（按 scaffold_id 注入）；L2 只表达结构 kind。


def _generic_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    """无保留 scaffold 时的通用环 kind(统一收敛为 alkane, 环系身份由 scaffold_id/fused_tree 承载)。"""
    mol = info["mol"]
    atoms = set(skeleton.atom_ids)
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms):
        n_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if set(ring) <= atoms)  # 芳香稠环(未注册)kind 收敛 alkane, 骨架身份由 fused_tree + scaffold_id 承载(L5 fused_namer 组装稠合名)。
        if n_rings >= 2:
            return "alkane"
        return None
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in atoms)
    return "alkane" if all_carbon else None  # 纯烃环统一 kind='alkane'（正交化：环系由 scaffold_id 承载，不饱和度由 double_bond/double_bonds 字段承载，命名由 chain_engine 动态加 cyclo 前缀）。


def _resolved_ring_kind(scaffold, info: dict, skeleton: ParentSkeleton) -> str | None:
    """按保留 scaffold 解析环 kind（苯/未注册稠环均收敛为 alkane）。"""
    if scaffold and scaffold.id != "carbocycle":
        if scaffold.id == "benzene":
            return "alkane"  # 苯环（纯烃芳香单环）kind 收敛为 alkane，环系由 scaffold_id="benzene" 承载（对齐环烷烃正交化）。
        if scaffold.id in ("fused", "fused_hetero"):
            return "alkane"  # 未注册稠环：kind 正交化收敛 alkane，骨架身份由 fused_tree + scaffold_id 承载（L5 fused_namer 组装稠合名）。
        return scaffold.id
    return _generic_ring_kind(info, skeleton)


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int, scaffold) -> str | None:
    """决定环骨架母体的 kind（radical/正交化 FG 类/结构 kind）。"""
    if selection.group_class is FunctionalGroupClass.RADICAL:
        if scaffold is None:  # 苯基取代基保留名（P-22.2.4）由 L5 radical worker 的 benzene variant 表达，L2 只给统一 kind。
            return None  # 未知杂环 scaffold：L5 无 -yl 词干可拼，显式失败而非当开链烷基错名。
        return "radical"
    if scaffold is not None and selection.group_class in _CHAIN_FG:  # 环 + 主 FG → FG 类别 kind（正交化）：苯/饱和环/稠环/杂环一律收敛，命名由 L5 chain_engine 通用词干引擎拼接（苯等保留名经 variant 特殊，无 variant 走通用名）。
        kind = _chain_kind(selection.group_class, count)
        if kind is not None:
            return kind
    if scaffold is not None and selection.group_class is FunctionalGroupClass.ALDEHYDE:
        return "aldehyde"  # 环上外环 -CHO 可多个同作主官能团（-carbaldehyde / -dicarbaldehyde，P-66.6.1.1.3）；aldehyde 不在 _MULTI_FG（开链二醛仍不支撑），此处仅环骨架放行，避免改变开链表达。
    return _resolved_ring_kind(scaffold, info, skeleton)


def _ring_fields(selection: PrincipalGroupSelection, occurrences) -> dict:
    """构造环主基团的 anchor/特征原子字段。"""
    anchors = _anchors(occurrences)
    fields = _principal_fields(selection.group_class, anchors) if _anchor_fields(selection.group_class) else {}
    characteristic = sorted({i for o in occurrences for i in o.characteristic_atoms})
    return {**fields, "principal_characteristic_atoms": characteristic,
            "principal_attachment_atoms": anchors,
            "principal_group_class": selection.group_class.value}


def _ring_fact_fields(fields: dict, facts: PrincipalExpressionFacts) -> dict:
    """合并附着原子字段，单附着酸/酯/酰胺/腈加 ring_attach_idx。"""
    attachments = sorted(facts.attachment_atoms)
    extra = {"principal_attachment_atoms": attachments}
    if len(attachments) == 1 and facts.group_class in (
        FunctionalGroupClass.ACID, FunctionalGroupClass.ESTER,
        FunctionalGroupClass.AMIDE, FunctionalGroupClass.NITRILE,
        FunctionalGroupClass.ALDEHYDE, FunctionalGroupClass.ACYL,
    ):
        extra["ring_attach_idx"] = attachments[0]  # ACYL：exocyclic 酰基头（苯甲酰/furan-2-carbonyl）的环附着原子位次，供 L4 在 locant_calc 计算 -carbonyl/benzoyl 词形所需 locant。
    return {**fields, **extra}


def _scaffold_fields(info: dict, skeleton: ParentSkeleton, facts=None, scaffold=None) -> dict:
    """解析并写入 scaffold 身份与表达能力字段；稠环拆解独立于 scaffold 身份。"""

    from namepredict.layer2.ring_expression_policy import supports_ring_expression
    if scaffold is None:
        from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
        scaffold = resolve_ring_scaffold(info, skeleton)  # print(scaffold)
    fields: dict = {}
    if scaffold:
        supported = supports_ring_expression(scaffold, facts) if facts else False
        match = None  # 保留 fused 模板匹配映射：L4 固定编号（standard_path）据此把模板原子映射到分子原子。
        from namepredict.layer2.ring_scaffold import _match_with_map, get_spec, hydrogenated_atoms
        spec = get_spec(scaffold.id)
        if spec and (spec.numbering.standard_path or spec.n_rings >= 2):  # 多环模板同样取 match：氢化衍生物需据此比对 mancude 参考算加氢位（P-31.2.2）
            hit = _match_with_map(info, skeleton.atom_ids)
            match = hit[1] if hit and hit[0] == scaffold.id else None

        fields = {"scaffold_id": scaffold.id, "scaffold_identity": scaffold,
                  "scaffold_match": match,
                  "typed_ring_expression_supported": supported}  # print({"scaffold_id": scaffold.id, "scaffold_identity": scaffold, "scaffold_match": match, "typed_ring_expression_supported": supported})
        if match:  # 加氢原子集：编号完成后由 L4 换算为 hydro 前缀位次
            hydro = hydrogenated_atoms(info["mol"], scaffold.id, match)
            if hydro:
                fields["hydro_atoms"] = hydro
    system = next((s for s in info.get("ring_systems") or []
                   if (s.get("atom_ids") or []) == list(skeleton.atom_ids)), None)  # 多环骨架附加稠环拆解结构（fused_tree 为 FusedNode 对象供 L5 稠合名组装）；拆解独立于 scaffold 身份：未注册系统 scaffold 解析为 None 时仍产出拆解树，供 L5 fused_namer 组装稠合名（P-25.3.2）。

    if system is not None and len(system.get("sssr_indices") or ()) >= 2:
        from namepredict.layer2.fused_system import decompose_fused_system
        node = decompose_fused_system(info, system)
        if node is not None:
            fields["fused_tree"] = node
    return fields


def ester_fields(info: dict, fields: dict) -> dict:
    """取首个酯的 o_idx 并写入酯字段（alkoxy_n 恒 0）。"""
    e = info["esters"][0]
    return {**fields, "o_idx": e["o_idx"], "alkoxy_n": 0}


def express_ring_principal(info: dict, selection: PrincipalGroupSelection,
                           skeleton: ParentSkeleton) -> dict | None:
    """环骨架：表达主基团并生成母体 dict（不支持返回 None）。"""
    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold  # 骨架原子集已确定：一次识别 scaffold，下游复用（不再重复调 _producer_id）。
    scaffold = resolve_ring_scaffold(info, skeleton)
    occurrences = _covered(selection, skeleton)
    kind = _ring_kind(info, selection, skeleton, len(occurrences), scaffold)
    if kind is None:
        return None
    facts = _facts(selection, skeleton, occurrences, info["mol"])
    fields = {**_ring_fact_fields(_ring_fields(selection, occurrences), facts),
              **_chain_unsat_fields(info, skeleton,
                                    _scaffold_fields(info, skeleton, facts, scaffold))}  # 补环内不饱和字段：kind 正交化后（醇/酮/纯烃环 → FG 类别/alkane），烯/炔由 double_bond(s)/triple_bond 字段承载（否则环烯酮/环烯醇/环烯烃烯丢失）。
    if facts.group_class is FunctionalGroupClass.ACID:
        fields = {**fields, **_expression_flags(selection, occurrences)}  # 环酸全阴离子补 anion 标志（链酸经 _chain_fields→_expression_flags 已设）；L5 据此转 -ate/-酸根，并让金属盐前缀（sodium …）能命中。
    if facts.group_class is FunctionalGroupClass.ESTER and facts.multiplicity == 1:
        fields = ester_fields(info, fields)
    if facts.group_class is FunctionalGroupClass.ACYL_HALIDE and facts.multiplicity == 1:
        fields = _chain_acyl_halide_fields(info, occurrences, fields)  # 环外酰卤（苯甲酰卤等）同样要卤素字段：hal_z 供 L5 选氟氯溴碘后缀，hal_idx 纳入母体原子。
    return _parent_dict(kind, skeleton, occurrences, fields, facts)


def _chain_fields(selection, occurrences) -> dict:
    """构造链主基团的 anchor 与表达标志字段。"""
    anchors = _anchors(occurrences)
    return {**_principal_fields(selection.group_class, anchors),
            **_expression_flags(selection, occurrences)}


def _unsat_bond_fields(dbs: list[dict], tbs: list[dict]) -> dict:
    """不饱和度 → 双键/三键字段字典（单数单键、列表多键，烯/炔可共存）。"""
    fields = {}
    if len(dbs) == 1:
        fields["double_bond"] = (dbs[0]["c1"], dbs[0]["c2"])
    elif len(dbs) >= 2:
        fields["double_bonds"] = [(d["c1"], d["c2"]) for d in dbs]
    if len(tbs) == 1:
        fields["triple_bond"] = (tbs[0]["c1"], tbs[0]["c2"])
    elif len(tbs) >= 2:
        fields["triple_bonds"] = [(t["c1"], t["c2"]) for t in tbs]
    return fields


def _chain_unsat_fields(info: dict, skeleton: ParentSkeleton, fields: dict) -> dict:
    """为骨架内的 C=C/C≡C 附加双键/三键位次字段。"""
    dbs, tbs = _chain_polys(info, set(skeleton.atom_ids))
    return {**fields, **_unsat_bond_fields(dbs, tbs)}


def _chain_phosphate_fields(info: dict, occurrences, fields: dict) -> dict | None:
    """L5 phosphate_names 消费的计数与盐元数据，盐门控不通过返回 None（同旧 build_phosphate_parent）：n_om>0 时若有碱金属须同数配对；中性酸/酯不允许带金属；完全无抗衡金属的游离磷酸根/磷酸酯阴离子放行。"""
    if len(occurrences) != 1:
        return None
    payload = occurrences[0].payload
    n_oh, n_om = int(payload.get("n_oh", 0)), int(payload.get("n_om", 0))
    salt = dict(info.get("salt") or {})
    if n_om > 0:
        if salt.get("metal") and int(salt.get("n_metal") or 0) != n_om:
            return None
    elif salt.get("metal"):
        return None
    return {**fields, "n_oh": n_oh, "n_om": n_om,
            "n_arms": int(payload.get("n_arms", 0)), "salt_meta": salt or None}


def _chain_ester_fields(info: dict, occurrences, fields: dict) -> dict:
    """L5 酯命名的酯烷氧基侧字段：o_idx 供 o_side 识别；仅严格线性给 alkoxy_n 保留名。"""
    if len(occurrences) != 1:
        return fields
    match = next((e for e in (info.get("esters") or [])
                  if e["c_idx"] in occurrences[0].characteristic_atoms), None)
    if match is None:
        return fields
    return {**fields, "o_idx": match["o_idx"], "alkoxy_n": 0}


# 表 2.1 单核母体氢化物（P-15.4.1）：杂原子锚点自由基的母体 free 名与元素。
# 去氢即标准取代基名（oxidane→hydroxy/oxy、azane→amino、sulfane→sulfanyl）。
_MONONUCLEAR_STEM: dict[int, tuple[str, str, str]] = {
    7: ("N", "azane", "氮烷"),
    8: ("O", "oxidane", "氧化烷"),
    16: ("S", "sulfane", "硫烷"),
}


def _mononuclear_radical(info: dict, skeleton: ParentSkeleton,
                         occurrences) -> tuple[ParentSkeleton, dict] | None:
    """杂原子锚点自由基收敛为单核氢化物母体骨架（表 2.1）：chain 限定单原子、注入 oxidane/azane 等 free 名，碳侧链留 L3 供 L5 经 free_to_yl 转标准取代基名；仅支持单一锚点。"""
    mol = info["mol"]
    anchors = sorted({i for o in occurrences for i in o.parent_anchors})
    if len(anchors) != 1:
        return None
    spec = _MONONUCLEAR_STEM.get(mol.GetAtomWithIdx(anchors[0]).GetAtomicNum())
    if spec is None:
        return None
    element, stem_en, stem_zh = spec
    new = replace(skeleton, atom_ids=(anchors[0],))
    return new, {"radical_anchor_element": element,
                 "stem_en": stem_en, "stem_zh": stem_zh}


def _chain_acyl_halide_fields(info: dict, occurrences, fields: dict) -> dict:
    """酰卤的卤素字段：hal_idx 供 parent_ownership 把卤素纳入母体原子（Cl 不作取代基）。"""
    if len(occurrences) != 1:
        return fields
    match = next((e for e in (info.get("acyl_chlorides") or [])
                  if e["c_idx"] in occurrences[0].characteristic_atoms), None)
    if match is None:
        return fields
    return {**fields, "hal_idx": match["hal_idx"], "hal_z": match["hal_z"]}


def express_chain_principal(info: dict, selection: PrincipalGroupSelection,
                            skeleton: ParentSkeleton) -> dict | None:
    """链骨架：表达主基团并生成母体 dict（不支持返回 None）。"""
    if skeleton.topology is not SkeletonTopology.ACYCLIC:
        return None

    occurrences = _covered(selection, skeleton)
    kind = _chain_kind(selection.group_class, len(occurrences))

    if kind is None:
        return None
    extra: dict = {}
    if kind == "radical":
        mono = _mononuclear_radical(info, skeleton, occurrences)
        if mono is not None:
            skeleton, extra = mono
    fields = _chain_unsat_fields(info, skeleton, {**_chain_fields(selection, occurrences), **extra})
    if kind == "ester":
        fields = _chain_ester_fields(info, occurrences, fields)
    elif kind == "acyl_halide":
        fields = _chain_acyl_halide_fields(info, occurrences, fields)
    elif kind == "phosphate":
        fields = _chain_phosphate_fields(info, occurrences, fields)
        if fields is None:
            return None
    return _parent_dict(kind, skeleton, occurrences, fields,
                        _facts(selection, skeleton, occurrences))


# ── 无主官能团（纯烃）表达：P-44.1 缺位时按拓扑分配 hydrocarbon kind ──

def _chain_polys(info: dict, atom_set: set[int]) -> tuple[list[dict], list[dict]]:
    """骨架内的 C=C / C≡C 条目（端点都在 atom_set 中）。"""
    dbs = [d for d in info.get("double_bonds") or [] if d["c1"] in atom_set and d["c2"] in atom_set]
    tbs = [t for t in info.get("triple_bonds") or [] if t["c1"] in atom_set and t["c2"] in atom_set]
    return dbs, tbs