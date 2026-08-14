# IUPAC: P-22.2.1 / P-25
# Layer: L2
"""Single source of truth: ScaffoldSpec drives kind_registry stems + retained ids."""
from __future__ import annotations

import pytest

from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.ring_scaffold import get_entry, registry
from namepredict.layer2.ring_scaffold import (
    FUSED56_SPECS,
    MONO_HETERO_SPECS,
    NAPH_FAMILY_SPECS,
    all_specs,
    fused56_kind_ids,
    get_spec,
    kind_ids_for,
    monohetero_kind_ids,
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
    for sid in (
        "anthracene", "quinazoline", "quinoxaline", "pyridine",
        "furan", "indole", "benzene",
    ):
        assert get_spec(sid) is not None, sid


def test_anthracene_meta() -> None:
    sp = get_spec("anthracene")
    assert sp is not None
    assert sp.n_rings == 3
    assert sp.ring == "carbo"
    assert sp.retained is True
    assert sp.stem_en == "anthracene" and sp.stem_zh == "蒽"


@pytest.mark.parametrize("sid", ["quinazoline", "quinoxaline"])
def test_benzodiazine_naph_family(sid: str) -> None:
    sp = get_spec(sid)
    assert sp is not None
    assert sp.naming_class == "benzodiazine"
    assert sp.n_rings == 2
    assert sp.ring == "hetero"
    assert len(sp.numbering.standard_path) == 10
    assert sp.numbering.mode == "naph_family"


def test_monohetero_table_covers_legacy() -> None:
    ids = {s.id for s in MONO_HETERO_SPECS}
    for sid in ("pyridine", "furan", "oxolane", "piperidine", "morpholine"):
        assert sid in ids
    assert kind_ids_for("monohetero") == monohetero_kind_ids()
    assert "pyridine" in monohetero_kind_ids()


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
    for sp in FUSED56_SPECS + NAPH_FAMILY_SPECS + MONO_HETERO_SPECS:
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
