"""Principal characteristic-group selection before P-44 parent selection."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass,
    FunctionalGroupInventory,
    FunctionalGroupOccurrence,
)


# Explicit precedence: first item is senior. Functional-prefix-only classes are absent.
_PRINCIPAL_ORDER = (
    FunctionalGroupClass.ACID,
    FunctionalGroupClass.SULFONIC_ACID,
    FunctionalGroupClass.ANHYDRIDE,
    FunctionalGroupClass.ESTER,
    FunctionalGroupClass.ACYL_HALIDE,
    FunctionalGroupClass.AMIDE,
    FunctionalGroupClass.NITRILE,
    FunctionalGroupClass.ALDEHYDE,
    FunctionalGroupClass.KETONE,
    FunctionalGroupClass.ALCOHOL,
    FunctionalGroupClass.THIOL,
    FunctionalGroupClass.AMINE,
)


@dataclass(frozen=True)
class PrincipalGroupSelection:
    group_class: FunctionalGroupClass
    occurrences: tuple[FunctionalGroupOccurrence, ...]


def select_principal_group(inventory: FunctionalGroupInventory) -> PrincipalGroupSelection | None:
    for group_class in _PRINCIPAL_ORDER:
        occurrences = inventory.occurrences(group_class)
        if occurrences:
            return PrincipalGroupSelection(group_class, occurrences)
    return None
