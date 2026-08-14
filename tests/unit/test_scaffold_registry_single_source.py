# IUPAC: P-22.2.1 / P-25
# Layer: L2
"""Single source of truth: ScaffoldSpec drives kind_registry stems + retained ids.

正交化后 _ALL_SPECS 只含 4 个保留 scaffold（benzene/pyridine/naphthalene/indole），
fused56/naph_family/monohetero 的固定位次标签已废弃（standard_path 为空）。
"""
from __future__ import annotations

import pytest

from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.ring_scaffold import (
    all_specs,
    fused56_kind_ids,
    get_entry,
    get_spec,
    kind_ids_for,
    monohetero_kind_ids,
    naph_kind_ids,
    registry,
)


# Specs that carry stems must match kind_registry.parent_names.
def _stem_specs():
    return [
        s for s in all_specs()
        if s.stem_en is not None and s.stem_zh is not None
    ]


@pytest.mark.parametrize("sp", _stem_specs(), ids=lambda s: s.id)
def test_spec_stem_matches_kind_registry(sp) -> None:
    assert kr.parent_names(sp.id) == (sp.stem_en, sp.stem_zh)


def test_required_scaffolds_in_spec() -> None:
    for sid in ("benzene", "pyridine", "naphthalene", "indole"):
        assert get_spec(sid) is not None, sid


@pytest.mark.parametrize("sid", ["benzene", "pyridine", "naphthalene", "indole"])
def test_spec_meta(sid: str) -> None:
    sp = get_spec(sid)
    assert sp is not None
    assert sp.retained is True
    assert sp.ring in ("carbo", "hetero")
    assert sp.n_rings in (1, 2)


def test_kind_ids_helpers_derive_from_specs() -> None:
    assert kind_ids_for("mono_carbo") == {"benzene"}
    assert kind_ids_for("monohetero") == {"pyridine"}
    assert kind_ids_for("naph_family") == {"naphthalene"}
    assert kind_ids_for("fused56") == {"indole"}
    assert fused56_kind_ids() == {"indole"}
    assert naph_kind_ids() == {"naphthalene"}
    assert monohetero_kind_ids() == {"pyridine"}


def test_fused56_ids_retained_in_kind_registry() -> None:
    for sid in fused56_kind_ids():
        m = kr.get(sid)
        assert m is not None, sid
        assert m.retained is True
        assert m.n_rings == 2


def test_monohetero_ids_retained_in_kind_registry() -> None:
    for sid in monohetero_kind_ids():
        m = kr.get(sid)
        assert m is not None, sid
        assert m.retained is True
        assert m.ring == "hetero"
        assert m.n_rings == 1


def test_no_bidirectional_drift_stem_specs() -> None:
    """Every retained Spec with stem is in kind_registry with same stems."""
    for sp in all_specs():
        if sp.stem_en is None or sp.stem_zh is None:
            continue
        assert kr.parent_names(sp.id) == (sp.stem_en, sp.stem_zh)


def test_retained_registry_ids_subset_of_spec() -> None:
    for sid in registry():
        assert get_spec(sid) is not None, sid
        entry = get_entry(sid)
        assert entry is not None
        sp = get_spec(sid)
        assert entry["en"] == sp.stem_en
        assert entry["zh"] == sp.stem_zh


def test_indole_parent_names_stable() -> None:
    assert kr.parent_names("indole") == ("1H-indole", "吲哚")
