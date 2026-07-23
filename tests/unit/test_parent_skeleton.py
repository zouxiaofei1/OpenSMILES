# IUPAC: P-44.1.1, P-44.3
# Layer: L1,L2
"""Topology-first principal skeleton enumeration."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import (
    SkeletonTopology, ParentSkeleton, enumerate_principal_skeletons,
    keep_max_principal_coverage, keep_max_ring_system_size, keep_p44_1_2, keep_p44_3, p44_3_key,
)
from namepredict.layer2.principal_selection import select_principal_group
from namepredict.layer2.principal_parent import select_principal_parent_skeletons
from namepredict.layer2.principal_expression import express_chain_principal


def _case(smiles: str, group_class: FunctionalGroupClass):
    info = analyze(preprocess(smiles))
    occurrences = info["fg_inventory"].occurrences(group_class)
    return occurrences, enumerate_principal_skeletons(info, occurrences)


def test_single_acid_finds_longest_related_open_chain() -> None:
    occurrences, result = _case("CC(C)CC(=O)O", FunctionalGroupClass.ACID)
    chains = [c for c in result.candidates if c.topology is SkeletonTopology.ACYCLIC]
    assert len(occurrences) == 1
    assert max(map(lambda c: len(c.atom_ids), chains)) == 4
    assert chains[0].covered_principal_ids == frozenset({"carboxyls:0"})


def test_two_acids_on_one_chain_are_both_covered() -> None:
    _, result = _case("O=C(O)CCC(=O)O", FunctionalGroupClass.ACID)
    assert max(map(lambda c: len(c.covered_principal_ids), result.candidates)) == 2
    assert not result.unsupported_ids


def test_three_branched_acids_produce_maximum_pair_candidates() -> None:
    _, result = _case("C(C(=O)O)(C(=O)O)C(=O)O", FunctionalGroupClass.ACID)
    counts = [len(c.covered_principal_ids) for c in result.candidates]
    assert max(counts) == 2
    assert len([n for n in counts if n == 2]) >= 2


def test_ring_exocyclic_acids_find_same_ring_system() -> None:
    _, result = _case("O=C(O)C1CCC(C(=O)O)CC1", FunctionalGroupClass.ACID)
    rings = [c for c in result.candidates if c.topology is SkeletonTopology.RING_SYSTEM]
    assert len(rings) == 1
    assert len(rings[0].atom_ids) == 6
    assert len(rings[0].covered_principal_ids) == 2


def test_principal_selection_ignores_functional_prefix_only_ether() -> None:
    info = analyze(preprocess("COc1ccccc1"))
    assert select_principal_group(info["fg_inventory"]) is None


def test_principal_selection_is_explicit_not_same_rank_order() -> None:
    info = analyze(preprocess("CC(=O)OCS(=O)(=O)O"))
    principal = select_principal_group(info["fg_inventory"])
    assert principal is not None
    assert principal.group_class is FunctionalGroupClass.SULFONIC_ACID


def test_p44_1_1_keeps_only_maximum_coverage() -> None:
    _, result = _case("C(C(=O)O)(C(=O)O)C(=O)O", FunctionalGroupClass.ACID)
    kept = keep_max_principal_coverage(result.candidates)
    assert kept
    assert {len(c.covered_principal_ids) for c in kept} == {2}


def test_p44_1_2_same_senior_atom_prefers_ring_over_chain() -> None:
    info = analyze(preprocess("CCCC1CCCCC1"))
    ring = ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(info["rings"][0]["atom_ids"]), frozenset())
    chain = ParentSkeleton(SkeletonTopology.ACYCLIC, (0, 1, 2, 3), frozenset())
    assert keep_p44_1_2(info["mol"], (chain, ring)) == (ring,)


def test_p44_1_2_does_not_decide_ring_vs_ring() -> None:
    info = analyze(preprocess("c1ccncc1-c2ccccc2"))
    rings = tuple(ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(r["atom_ids"]), frozenset()) for r in info["rings"])
    kept = keep_p44_1_2(info["mol"], rings)
    assert kept == rings


def test_largest_related_ring_system_wins_after_coverage() -> None:
    small = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (0, 1, 2), frozenset({"x"}))
    large = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (3, 4, 5, 6), frozenset({"x"}))
    assert keep_max_ring_system_size((small, large)) == (large,)


def test_equal_size_ring_systems_remain_tied() -> None:
    first = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (0, 1, 2), frozenset({"x"}))
    second = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (3, 4, 5), frozenset({"x"}))
    assert keep_max_ring_system_size((first, second)) == (first, second)


def test_principal_coverage_precedes_ring_system_size() -> None:
    small = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (0, 1, 2), frozenset({"x", "y"}))
    large = ParentSkeleton(SkeletonTopology.RING_SYSTEM, (3, 4, 5, 6), frozenset({"x"}))
    covered = keep_max_principal_coverage((small, large))
    assert keep_max_ring_system_size(covered) == (small,)


def test_p44_3_longer_carbon_chain_precedes_unsaturation() -> None:
    info = analyze(preprocess("CCCC.C=C"))
    long = ParentSkeleton(SkeletonTopology.ACYCLIC, (0, 1, 2, 3), frozenset())
    short = ParentSkeleton(SkeletonTopology.ACYCLIC, (4, 5), frozenset())
    assert p44_3_key(info["mol"], long) > p44_3_key(info["mol"], short)
    assert keep_p44_3(info["mol"], (short, long)) == (long,)


def test_rule_driven_entry_never_calls_fg_producers(monkeypatch) -> None:
    from namepredict.layer2 import kind_registry

    monkeypatch.setattr(kind_registry, "fg_try_fns", lambda: (_ for _ in ()).throw(AssertionError("blind try")))
    selected = select_principal_parent_skeletons(analyze(preprocess("O=C(O)CCC(=O)O")))
    assert selected.principal is not None
    assert len(selected.skeletons.candidates[0].covered_principal_ids) == 2


def test_expression_uses_selected_skeleton_without_producer() -> None:
    selected = select_principal_parent_skeletons(analyze(preprocess("O=C(O)CCC(=O)O")))
    parent = express_chain_principal(selected.principal, selected.skeletons.candidates[0])
    assert parent["kind"] == "diacid"
    assert parent["principal_group_count"] == 2
    assert parent["covered_principal_ids"] == ("carboxyls:0", "carboxyls:1")


def test_production_path_has_no_legacy_fg_candidates() -> None:
    from namepredict.layer2 import candidates

    assert not hasattr(candidates, "_fg_candidates")
    parents = candidates._collect_candidates(analyze(preprocess("O=C(O)CCC(=O)O")))
    assert parents[0]["kind"] == "diacid"
    assert parents[0]["principal_group_count"] == 2


def test_no_principal_group_defers_to_hydrocarbon_selection() -> None:
    selected = select_principal_parent_skeletons(analyze(preprocess("COc1ccccc1")))
    assert selected.principal is None
    assert selected.skeletons is None


def test_remote_ring_is_not_related_to_acid() -> None:
    _, result = _case("O=C(O)CCc1ccccc1", FunctionalGroupClass.ACID)
    rings = [c for c in result.candidates if c.topology is SkeletonTopology.RING_SYSTEM]
    assert not rings
