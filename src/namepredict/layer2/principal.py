"""P-41 类 / P-43 表达式元数据 + 主基团选择。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping

from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass as FG,
    FunctionalGroupInventory,
    FunctionalGroupOccurrence,
)


@dataclass(frozen=True, order=True)
class PrincipalPriority:
    p41_class: int
    p43_path: tuple[int, ...] = ()


class PrincipalExpression(str, Enum):
    SUFFIX = "suffix"
    PREFIX_ONLY = "prefix_only"
    LEGACY_COMPAT = "legacy_compat"


@dataclass(frozen=True)
class PrincipalFeatureSpec:
    priority: PrincipalPriority
    expression: PrincipalExpression
    compatibility_rank: int = 0
    anchor_fields: tuple[str, str] | None = None


def _suffix(p41_class: int, rank: int, *path: int,
            anchor_fields: tuple[str, str] | None = None) -> PrincipalFeatureSpec:
    """构造 SUFFIX 型主官能团规格。"""
    priority = PrincipalPriority(p41_class, path)
    return PrincipalFeatureSpec(priority, PrincipalExpression.SUFFIX, rank, anchor_fields)


PRINCIPAL_REGISTRY: dict[FG, PrincipalFeatureSpec] = {
    FG.RADICAL: _suffix(1, 1, anchor_fields=("radical_c_idx", "radical_c_idxs")),
    FG.ACID: _suffix(7, 14, 1, anchor_fields=("cooh_c_idx", "cooh_c_idxs")),
    FG.ANHYDRIDE: _suffix(8, 12),
    FG.ESTER: _suffix(9, 11, anchor_fields=("ester_c_idx", "ester_c_idxs")),
    FG.ACYL_HALIDE: _suffix(10, 10),
    FG.AMIDE: _suffix(11, 9, anchor_fields=("amide_c_idx", "amide_c_idxs")),
    FG.NITRILE: _suffix(14, 8, anchor_fields=("nitrile_c_idx", "nitrile_c_idxs")),
    FG.ALDEHYDE: _suffix(15, 7, anchor_fields=("aldehyde_c_idx", "aldehyde_c_idxs")),
    FG.KETONE: _suffix(16, 6, anchor_fields=("ketone_c_idx", "ketone_c_idxs")),
    FG.ALCOHOL: _suffix(17, 5, 1, anchor_fields=("oh_c_idx", "oh_c_idxs")),
    FG.THIOL: _suffix(17, 4, 2),
    FG.AMINE: _suffix(19, 3, anchor_fields=("amine_c_idx", "amine_c_idxs")),
    FG.ISOCYANATE: PrincipalFeatureSpec(PrincipalPriority(41), PrincipalExpression.LEGACY_COMPAT, 8),
    FG.ISOTHIOCYANATE: PrincipalFeatureSpec(PrincipalPriority(41), PrincipalExpression.LEGACY_COMPAT, 8),
    FG.SULFIDE: PrincipalFeatureSpec(PrincipalPriority(41, (2,)), PrincipalExpression.LEGACY_COMPAT, 2),
    FG.ETHER: PrincipalFeatureSpec(PrincipalPriority(41, (1,)), PrincipalExpression.PREFIX_ONLY),
    # Legacy_compat：仅取 kind_registry fg_rank 投影的 compatibility_rank，永不作主基团（principal_spec() 以 SUFFIX 为门槛）。
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
