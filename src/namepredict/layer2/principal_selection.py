"""Principal characteristic-group selection before P-44 parent selection."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass,
    FunctionalGroupInventory,
    FunctionalGroupOccurrence,
)
from namepredict.layer2.principal_registry import (
    PRINCIPAL_REGISTRY,
    PrincipalFeatureSpec,
    principal_spec,
)


@dataclass(frozen=True)
class PrincipalGroupSelection:
    group_class: FunctionalGroupClass
    occurrences: tuple[FunctionalGroupOccurrence, ...]


def select_principal_group(
    inventory: FunctionalGroupInventory,
    registry: Mapping[FunctionalGroupClass, PrincipalFeatureSpec] = PRINCIPAL_REGISTRY,
) -> PrincipalGroupSelection | None:
    classes = (entry.group_class for entry in inventory.entries)
    eligible = (group_class for group_class in set(classes) if principal_spec(group_class, registry))
    selected = min(eligible, key=lambda group_class: principal_spec(group_class, registry).priority, default=None)
    occurrences = inventory.occurrences(selected) if selected else ()
    return PrincipalGroupSelection(selected, occurrences) if selected else None
