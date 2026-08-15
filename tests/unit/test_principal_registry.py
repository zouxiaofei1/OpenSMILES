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
from namepredict.layer2.parent_candidate import _kind_rank
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
    custom[FG.ISOCYANATE] = PrincipalFeatureSpec(
        PrincipalPriority(7, (3, 2)), PrincipalExpression.SUFFIX,
    )
    selected = select_principal_group(_inventory(FG.ISOCYANATE, FG.AMIDE), custom)
    assert selected is not None and selected.group_class == FG.ISOCYANATE


@pytest.mark.parametrize("kind,rank", [
    ("sulfide", 2), ("isocyanate", 8),
    # Migrated from _CHAIN_FG fallback to PRINCIPAL_REGISTRY projection:
    # ranks must be preserved exactly.
    ("tetraalkylammonium", 0),
])
def test_kind_rank_projection_preserved(kind: str, rank: int) -> None:
    assert _kind_rank(kind) == rank


@pytest.mark.parametrize("kind", ["carboxylate", "formamide_like", "amine_oxide", "unknown_one"])
def test_kind_rank_never_uses_substring_guessing(kind: str) -> None:
    assert _kind_rank(kind) == 0


def test_chain_fg_anchor_fields_present() -> None:
    """链式表达能产生 kind 的 FG 都必须有 payload 锚点字段（_FIELDS 并入主表后的守卫）。"""
    from namepredict.layer2.principal_expression import _CHAIN_FG, _anchor_fields
    assert all(_anchor_fields(gc) is not None for gc in _CHAIN_FG)
