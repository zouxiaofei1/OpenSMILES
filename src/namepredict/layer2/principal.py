"""P-41 类 / P-43 表达式元数据 + 主基团选择。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping

from namepredict.layer1.fg_registry import FgSpec, FG_SPECS, OXO_ACID_P41, oxoacid_is_acid
from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass as FG,
    FunctionalGroupInventory,
    FunctionalGroupOccurrence,
)


@dataclass(frozen=True, order=True)
class PrincipalPriority:
    """主基团优先级：P-41 类号 + P-43 路径（可比较排序）。"""
    p41_class: int
    p43_path: tuple[int, ...] = ()


class PrincipalExpression(str, Enum):
    """主基团的表达方式（作后缀）。"""
    SUFFIX = "suffix"


@dataclass(frozen=True)
class PrincipalFeatureSpec:
    """主官能团规格：优先级、表达方式与锚点字段名。"""
    priority: PrincipalPriority
    expression: PrincipalExpression
    anchor_fields: tuple[str, str] | None = None


def _spec_from_fg(sp: FgSpec) -> PrincipalFeatureSpec:
    """由 FG 注册表条目构造主官能团规格（优先级/表达/锚点字段）。"""
    return PrincipalFeatureSpec(
        PrincipalPriority(sp.p41, sp.path),
        PrincipalExpression(sp.expr),
        sp.parent_anchor_fields,
    )

PRINCIPAL_REGISTRY: dict[FG, PrincipalFeatureSpec] = { 
    FG(sp.fg): _spec_from_fg(sp) for sp in FG_SPECS if sp.p41
}


def feature_spec(group_class: FG, registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY) -> PrincipalFeatureSpec | None:
    """查 registry 返回基团类对应的规格（无则 None）。"""
    return registry.get(group_class)

@dataclass(frozen=True)
class PrincipalGroupSelection:
    """选中的主官能团类及其全部 occurrence。"""
    group_class: FG
    occurrences: tuple[FunctionalGroupOccurrence, ...]


_CATION_ANION_GATED = 99  # 阴离子在场时阳离子让位（表 4.1 类 4 > 类 6）：置底即不再作母体
_ANION_OS_P41 = 4  # 表 4.1 类 4 阴离子：氧/硫负离子作特征基团时高于酸（类 7）


def _has_charged_acid(inventory: FunctionalGroupInventory, mol) -> bool:
    """是否有阴离子型酸：酸根仍作母体，中性酸不阻断。"""
    if mol is None:
        return False
    for occ in inventory.occurrences(FG.ACID):  # 羧酸根：阴离子氧在 surr_idx 中
        if any(mol.GetAtomWithIdx(int(i)).GetFormalCharge() < 0
               for i in (occ.payload.get("surr_idx") or ())):
            return True
    for occ in inventory.occurrences(FG.OXOACID):  # 含氧酸根：n_om 即去质子臂数
        if int(occ.payload.get("n_om") or 0) > 0:
            return True
    return False


def _effective_priority(group_class: FG, spec: PrincipalFeatureSpec,
                        inventory: FunctionalGroupInventory, mol=None) -> PrincipalPriority:
    """候选类的实际 P-41 优先级：含氧酸式时升到类 7。"""
    if group_class in (FG.ALCOHOL, FG.THIOL):
        occurrences = inventory.occurrences(group_class)
        if (occurrences and all(o.payload.get("anion_os") for o in occurrences)
                and not _has_charged_acid(inventory, mol)):  # 酸根在场时仍由酸作母体（数据实测）
            return PrincipalPriority(_ANION_OS_P41, spec.priority.p43_path)  # 酚盐/醇盐/硫醇盐
    if group_class is FG.CATION and inventory.has_anion:
        return PrincipalPriority(_CATION_ANION_GATED, spec.priority.p43_path)
    if group_class is not FG.OXOACID:
        return spec.priority
    occurrences = inventory.occurrences(group_class)
    if occurrences and all(oxoacid_is_acid(o.payload) for o in occurrences):
        return PrincipalPriority(OXO_ACID_P41, spec.priority.p43_path)
    return spec.priority


def select_principal_group(
    inventory: FunctionalGroupInventory,
    registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY,
    mol=None,
) -> PrincipalGroupSelection | None:
    """按优先级选主官能团类并取全部 occurrence。"""
    classes = (entry.group_class for entry in inventory.entries if not entry.demoted)  # 降级叶（carboxy/cyano）不再作主基团候选
    eligible = (group_class for group_class in set(classes) if feature_spec(group_class, registry))
    selected = min(eligible, key=lambda group_class: _effective_priority(
        group_class, feature_spec(group_class, registry), inventory, mol), default=None)
    occurrences = inventory.occurrences(selected) if selected else ()
    return PrincipalGroupSelection(selected, occurrences) if selected else PrincipalGroupSelection(FG.NONE, occurrences)
