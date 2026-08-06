"""P-41 class / P-43 expression metadata + principal group selection."""
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


def _suffix(p41_class: int, rank: int, *path: int) -> PrincipalFeatureSpec:
    priority = PrincipalPriority(p41_class, path)
    return PrincipalFeatureSpec(priority, PrincipalExpression.SUFFIX, rank)


PRINCIPAL_REGISTRY: dict[FG, PrincipalFeatureSpec] = {
    FG.ACID: _suffix(7, 14, 1),
    FG.SULFONIC_ACID: _suffix(7, 13, 4),
    FG.ANHYDRIDE: _suffix(8, 12),
    FG.ESTER: _suffix(9, 11),
    FG.ACYL_HALIDE: _suffix(10, 10),
    FG.AMIDE: _suffix(11, 9),
    FG.NITRILE: _suffix(14, 8),
    FG.ALDEHYDE: _suffix(15, 7),
    FG.KETONE: _suffix(16, 6),
    FG.ALCOHOL: _suffix(17, 5, 1),
    FG.THIOL: _suffix(17, 4, 2),
    FG.AMINE: _suffix(19, 3),
    FG.CARBAMATE: PrincipalFeatureSpec(PrincipalPriority(9), PrincipalExpression.LEGACY_COMPAT, 11),
    FG.HYDRAZINE: PrincipalFeatureSpec(PrincipalPriority(21), PrincipalExpression.LEGACY_COMPAT, 4),
    FG.ISOCYANATE: PrincipalFeatureSpec(PrincipalPriority(41), PrincipalExpression.LEGACY_COMPAT, 8),
    FG.ISOTHIOCYANATE: PrincipalFeatureSpec(PrincipalPriority(41), PrincipalExpression.LEGACY_COMPAT, 8),
    FG.SULFIDE: PrincipalFeatureSpec(PrincipalPriority(41, (2,)), PrincipalExpression.LEGACY_COMPAT, 2),
    FG.SULFOXIDE: PrincipalFeatureSpec(PrincipalPriority(41, (3,)), PrincipalExpression.LEGACY_COMPAT, 2),
    FG.SULFONE: PrincipalFeatureSpec(PrincipalPriority(41, (4,)), PrincipalExpression.LEGACY_COMPAT, 6),
    FG.ETHER: PrincipalFeatureSpec(PrincipalPriority(41, (1,)), PrincipalExpression.PREFIX_ONLY),
}


def feature_spec(group_class: FG, registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY) -> PrincipalFeatureSpec | None:
    return registry.get(group_class)


def principal_spec(group_class: FG, registry: Mapping[FG, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY) -> PrincipalFeatureSpec | None:
    spec = feature_spec(group_class, registry)
    return spec if spec and spec.expression is PrincipalExpression.SUFFIX else None


def legacy_rank(group_class: FG | None) -> int:
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
    classes = (entry.group_class for entry in inventory.entries)
    eligible = (group_class for group_class in set(classes) if principal_spec(group_class, registry))
    selected = min(eligible, key=lambda group_class: principal_spec(group_class, registry).priority, default=None)
    occurrences = inventory.occurrences(selected) if selected else ()
    return PrincipalGroupSelection(selected, occurrences) if selected else None
