# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4
"""Retained fused aza scaffolds (indole/bim/quinoline/naph) as ScaffoldSpec."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.benzimidazole import _try_benzimidazole_parent
from namepredict.layer2.indole import _try_indole_parent
from namepredict.layer2.quinoline import _try_quinoline_parent
from namepredict.layer2.scaffold.specs import (
    FUSED56_LABELS,
    NAPH_LABELS,
    ScaffoldSpec,
    fused56_kind_ids,
    get_spec,
    kind_ids_for,
    naph_kind_ids,
)
from namepredict.layer4.locants.adapt import (
    FUSED56_KINDS,
    INDOLE_LABELS,
    INDOLE_ORIENT_KINDS,
    NAPH_KINDS,
    NAPH_LABELS as L4_NAPH_LABELS,
    Q_KINDS,
    plan_from_chain,
)
from namepredict.namer import SMILESNNamer

_F56_LABELS = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")
_NAPH_LABELS = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")

# Positive retained parents: kind, smiles, en, zh
POS_CASES = [
    ("indole", "c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("benzimidazole", "c1ccc2[nH]cnc2c1", "1H-benzimidazole", "1H-苯并咪唑"),
    ("quinoline", "c1ccc2ncccc2c1", "quinoline", "喹啉"),
    ("indazole", "c1ccc2[nH]ncc2c1", "1H-indazole", "1H-吲唑"),
    ("naphthalene", "c1ccc2ccccc2c1", "naphthalene", "萘"),
]

# Negative: wrong naming_class / not these families
NEG_CASES = [
    ("benzene", None),  # no retained fused Spec
    ("furan", None),
    ("cycloalkane", "carbocycle"),  # carbo free, not fused56/naph
]


def _mol(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    return mol


@pytest.mark.parametrize(
    "sid,nclass,n_labs",
    [
        ("indole", "fused56", 9),
        ("indolecarboxylic", "fused56", 9),
        ("indazole", "fused56", 9),
        ("indazolecarbonitrile", "fused56", 9),
        ("indazolecarbaldehyde", "fused56", 9),
        ("benzimidazole", "fused56", 9),
        ("benzimidazolamine", "fused56", 9),
        ("quinoline", "naph_family", 10),
        ("isoquinoline", "naph_family", 10),
        ("quinolinol", "naph_family", 10),
        ("quinolinecarboxylic", "naph_family", 10),
        ("naphthalene", "naph_family", 10),
        ("naphthalenecarboxylic", "naph_family", 10),
    ],
)
def test_retained_fused_specs_registered(sid: str, nclass: str, n_labs: int) -> None:
    sp = get_spec(sid)
    assert isinstance(sp, ScaffoldSpec)
    assert sp.id == sid
    assert sp.naming_class == nclass
    assert sp.n_rings == 2
    assert sp.retained is True
    assert len(sp.numbering.standard_path) == n_labs
    if nclass == "fused56":
        assert sp.numbering.mode == "fused56_fixed"
        assert sp.numbering.standard_path == _F56_LABELS
        assert sp.ring == "hetero"
    else:
        assert sp.numbering.mode in ("naph_family", "naph_fixed")
        assert sp.numbering.standard_path == _NAPH_LABELS


def test_labels_constants() -> None:
    assert FUSED56_LABELS == _F56_LABELS == INDOLE_LABELS
    assert NAPH_LABELS == _NAPH_LABELS == L4_NAPH_LABELS
    assert len(NAPH_LABELS) == 10
    assert "4a" in NAPH_LABELS and "8a" in NAPH_LABELS


def test_kind_id_helpers() -> None:
    f56 = fused56_kind_ids()
    assert "indole" in f56 and "benzimidazole" in f56
    assert "benzofuran" in f56  # prior fused56 retained
    assert "quinoline" not in f56
    naph = naph_kind_ids()
    assert "quinoline" in naph and "naphthalene" in naph
    assert "indole" not in naph
    assert kind_ids_for("fused56") == f56
    assert kind_ids_for("naph_family") == naph


@pytest.mark.parametrize("kind,nclass", NEG_CASES)
def test_negative_not_wrong_class(kind: str, nclass: str | None) -> None:
    sp = get_spec(kind)
    if nclass is None:
        assert sp is None or sp.naming_class not in ("fused56", "naph_family")
        return
    assert sp is not None
    assert sp.naming_class == nclass
    assert sp.naming_class not in ("fused56", "naph_family") or kind == "never"


def test_plan_from_chain_indole_9() -> None:
    plan = plan_from_chain(list(range(9)), "indole")
    assert plan is not None
    assert plan.labels == _F56_LABELS
    assert plan.scaffold_id == "indole"


def test_plan_from_chain_quinoline_10() -> None:
    plan = plan_from_chain(list(range(10)), "quinoline")
    assert plan is not None
    assert plan.labels == _NAPH_LABELS
    assert plan.scaffold_id == "quinoline"


def test_plan_from_chain_naphthalene_10() -> None:
    plan = plan_from_chain(list(range(10)), "naphthalene")
    assert plan is not None
    assert plan.labels == _NAPH_LABELS


def test_parent_indole_scaffold_id() -> None:
    parent = _try_indole_parent(analyze(_mol("c1ccc2[nH]ccc2c1")))
    assert parent is not None
    assert parent.get("kind") == "indole"
    assert parent.get("scaffold_id") == "indole"
    assert len(parent.get("chain") or []) == 9


def test_parent_bim_scaffold_id() -> None:
    parent = _try_benzimidazole_parent(analyze(_mol("c1ccc2[nH]cnc2c1")))
    assert parent is not None
    assert parent.get("scaffold_id") == "benzimidazole"


def test_parent_quinoline_scaffold_id() -> None:
    parent = _try_quinoline_parent(analyze(_mol("c1ccc2ncccc2c1")))
    assert parent is not None
    assert parent.get("scaffold_id") == "quinoline"
    assert len(parent.get("chain") or []) == 10


@pytest.mark.parametrize("smiles,kind,en,zh", [
    (s, k, e, z) for k, s, e, z in POS_CASES
])
def test_e2e_retained_names_stable(smiles, kind, en, zh) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_l2_l4_fused56_kind_contract() -> None:
    """L2 Spec fused56 ids must equal L4 FUSED56_KINDS (no L4→L2 import)."""
    assert fused56_kind_ids() == FUSED56_KINDS
    assert "indole" in FUSED56_KINDS
    assert "quinoline" not in FUSED56_KINDS


def test_l2_l4_naph_kind_contract() -> None:
    """L2 Spec naph_family ids must equal L4 NAPH_KINDS ∪ Q_KINDS."""
    assert naph_kind_ids() == (NAPH_KINDS | Q_KINDS)
    assert "quinoline" in Q_KINDS
    assert "naphthalene" in NAPH_KINDS


def test_l4_orient_kinds_union() -> None:
    """INDOLE_ORIENT_KINDS = fused56 + quinoline family (fixed chain orient)."""
    expect = frozenset(FUSED56_KINDS | Q_KINDS)
    assert frozenset(INDOLE_ORIENT_KINDS) == expect
    assert "naphthalene" not in expect  # naph uses dedicated orienter


def test_l2_l4_labels_contract() -> None:
    assert FUSED56_LABELS == INDOLE_LABELS
    assert NAPH_LABELS == L4_NAPH_LABELS
