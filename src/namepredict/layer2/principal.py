"""P-41 类 / P-43 表达式元数据 + 主基团选择。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping

from namepredict.layer1.fg_registry import FgSpec, FG_SPECS
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
    """主基团的表达方式（作后缀/仅作前缀/仅兼容等级）。"""
    SUFFIX = "suffix"
    PREFIX_ONLY = "prefix_only"
    LEGACY_COMPAT = "legacy_compat"


@dataclass(frozen=True)
class PrincipalFeatureSpec:
    """主官能团规格：优先级、表达方式、兼容等级与锚点字段名。"""
    priority: PrincipalPriority
    expression: PrincipalExpression
    compatibility_rank: int = 0
    anchor_fields: tuple[str, str] | None = None


def _spec_from_fg(sp: FgSpec) -> PrincipalFeatureSpec:
    """由 FG 注册表条目构造主官能团规格（优先级/表达/兼容等级/锚点字段）。"""
    return PrincipalFeatureSpec(
        PrincipalPriority(sp.p41, sp.path),
        PrincipalExpression(sp.expr),
        sp.compat,
        sp.parent_anchor_fields,
    )


# 主官能团规格唯一事实来源在 fg_registry.FG_SPECS；此处派生。
# Legacy_compat：仅取 compatibility_rank（legacy_rank），永不作主基团（principal_spec() 以 SUFFIX 为门槛）。
PRINCIPAL_REGISTRY: dict[FG, PrincipalFeatureSpec] = {
    FG(sp.fg): _spec_from_fg(sp) for sp in FG_SPECS if sp.p41
}


def feature_spec(group_class: FG, registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY) -> PrincipalFeatureSpec | None:
    """查 registry 返回基团类对应的规格（无则 None）。"""
    return registry.get(group_class)


def principal_spec(group_class: FG, registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY) -> PrincipalFeatureSpec | None:
    """仅返回 SUFFIX 型规格（Legacy 不作主基团）。"""
    spec = feature_spec(group_class, registry)
    return spec if spec and spec.expression is PrincipalExpression.SUFFIX else None


def legacy_rank(group_class: FG | None) -> int:
    """取基团类的兼容等级（无规格为 0）。"""
    spec = feature_spec(group_class) if group_class else None
    return spec.compatibility_rank if spec else 0


@dataclass(frozen=True)
class PrincipalGroupSelection:
    """选中的主官能团类及其全部 occurrence。"""
    group_class: FG
    occurrences: tuple[FunctionalGroupOccurrence, ...]


def select_principal_group(
    inventory: FunctionalGroupInventory,
    registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY,
) -> PrincipalGroupSelection | None:
    """按优先级选主官能团类并取全部 occurrence。"""
    classes = (entry.group_class for entry in inventory.entries)
    eligible = (group_class for group_class in set(classes) if principal_spec(group_class, registry))
    selected = min(eligible, key=lambda group_class: principal_spec(group_class, registry).priority, default=None)
    occurrences = inventory.occurrences(selected) if selected else ()
    return PrincipalGroupSelection(selected, occurrences) if selected else PrincipalGroupSelection(FG.NONE, occurrences)
