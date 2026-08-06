# IUPAC: P-41/P-43
# Layer: L2
"""Principal characteristic-group metadata is the sole priority authority."""
from __future__ import annotations

import pytest

from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass as FG,
    FunctionalGroupInventory,
    FunctionalGroupOccurrence,
)
from namepredict.layer2 import kind_registry
from namepredict.layer2.principal import (
    PRINCIPAL_REGISTRY,
    PrincipalExpression,
    PrincipalFeatureSpec,
    PrincipalPriority,
    select_principal_group,
)

# Representative structures for the P-41/P-43 boundaries under test.
CASES = [
    ("CC(=O)OCC(=O)O", FG.ACID, FG.ESTER),
    ("CS(=O)(=O)OC(=O)OC", FG.SULFONIC_ACID, FG.ANHYDRIDE),
    ("N#CCC=O", FG.NITRILE, FG.ALDEHYDE),
    ("COC", None, FG.ETHER),
]


def _occ(group_class: FG) -> FunctionalGroupOccurrence:
    return FunctionalGroupOccurrence(group_class.value, group_class, frozenset(), frozenset(), {})


def _inventory(*classes: FG) -> FunctionalGroupInventory:
    return FunctionalGroupInventory(tuple(_occ(group_class) for group_class in classes))


@pytest.mark.parametrize("smiles,expected,other", CASES)
def test_principal_boundaries(smiles: str, expected: FG | None, other: FG) -> None:
    del smiles
    selected = select_principal_group(_inventory(*(() if expected is None else (expected,)), other))
    assert (selected.group_class if selected else None) == expected


def test_registry_supports_new_p41_class_and_p43_subpath() -> None:
    custom = dict(PRINCIPAL_REGISTRY)
    custom[FG.PHOSPHONIC] = PrincipalFeatureSpec(
        PrincipalPriority(7, (3, 2)), PrincipalExpression.SUFFIX,
    )
    selected = select_principal_group(_inventory(FG.PHOSPHONIC, FG.AMIDE), custom)
    assert selected is not None and selected.group_class == FG.PHOSPHONIC


def test_legacy_kind_rank_is_registry_projection() -> None:
    assert kind_registry.fg_rank("acid") > kind_registry.fg_rank("ester")
    assert kind_registry.fg_rank("sulfonic_acid") > kind_registry.fg_rank("anhydride")
    assert kind_registry.fg_rank("nitrile") > kind_registry.fg_rank("aldehyde")
    assert kind_registry.fg_rank("ether") == 0
    assert kind_registry.has_principal_fg("ether") == 0


@pytest.mark.parametrize("kind,rank", [
    ("carbamate", 11), ("hydrazine", 4), ("sulfone", 6),
    ("sulfide", 2), ("sulfoxide", 2), ("isocyanate", 8),
    # Migrated from _CHAIN_FG fallback to PRINCIPAL_REGISTRY projection:
    # ranks must be preserved exactly.
    ("carbonate", 11), ("sulfonate", 11), ("sulfonyl_chloride", 10),
    ("urea", 9), ("guanidine", 9), ("sulfonamide", 9),
    ("tetraalkylammonium", 0), ("phosphate", 2), ("phosphonic", 2),
])
def test_unmigrated_legacy_ranks_are_preserved(kind: str, rank: int) -> None:
    assert kind_registry.fg_rank(kind) == rank


@pytest.mark.parametrize("kind", ["carboxylate", "formamide_like", "amine_oxide", "unknown_one"])
def test_kind_projection_never_uses_substring_guessing(kind: str) -> None:
    assert kind_registry.fg_rank(kind) == 0
