# IUPAC: P-44 / architecture
# Layer: L2,L5
"""ParentKind registry: single source for scoring meta + L5 parent stems.

正交化后 kind 只表达 FG 类别 / 保留 scaffold；环系维度由 scaffold_id 承载，
命名 kind（cycloalcohol/benzoic/…）由 L5 typed_kinds 决定，不注册在此。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.scoring import _score_parent
from namepredict.namer import SMILESNNamer

# 保留 scaffold kind（ring_scaffold._ALL_SPECS 4 个）必须有 stem + ring meta。
_RETAINED_KINDS = {
    "benzene": ("benzene", "苯"),
    "pyridine": ("pyridine", "吡啶"),
    "naphthalene": ("naphthalene", "萘"),
    "indole": ("1H-indole", "吲哚"),
}

@pytest.mark.parametrize("kind,names", _RETAINED_KINDS.items())
def test_retained_kind_stems(kind: str, names: tuple[str, str]) -> None:
    m = kr.get(kind)
    assert m is not None
    assert kr.parent_names(kind) == names
    assert m.retained is True


def test_ether_not_registered() -> None:
    # ether 非主官能团（compat 0），不作为 scaffold 词干注册。
    assert kr.get("ether") is None


def test_unknown_kind_negative() -> None:
    assert kr.get("no_such") is None
    assert kr.parent_names("no_such") is None
    assert kr.is_hetero_ring("no_such") == 0
    assert kr.is_carbo_ring("no_such") == 0
    assert kr.n_rings_of("no_such") == 0
    assert kr.retained_bonus("no_such") == 0


# 正交化契约：环 + 主 FG → FG 类别 kind + scaffold_id；无主 FG 保留 scaffold 结构 kind。
_ORTHO_CASES = [
    ("Oc1ccccc1", "alcohol", "benzene"),
    ("Nc1ccccc1", "amine", "benzene"),
    ("O=C(O)c1ccccc1", "acid", "benzene"),
    ("n1ccccc1", "pyridine", "pyridine"),
    ("c1ccc2ccccc2c1", "naphthalene", "naphthalene"),
    ("c1ccc2[nH]ccc2c1", "indole", "indole"),
]


@pytest.mark.parametrize("smiles,kind,scaffold", _ORTHO_CASES)
def test_select_parent_orthogonalized_kind(smiles: str, kind: str, scaffold: str) -> None:
    mol = preprocess(smiles)
    assert mol is not None
    parent = select_parent(analyze(mol))
    assert parent["kind"] == kind
    assert parent.get("scaffold_id") == scaffold


def test_select_parent_preloads_registry_stem() -> None:
    """有 spec stem 的无 FG 环 scaffold：select_parent 预填 kind_registry 词干。"""
    for smiles, kind in (("n1ccccc1", "pyridine"), ("c1ccc2ccccc2c1", "naphthalene")):
        mol = preprocess(smiles)
        assert mol is not None
        parent = select_parent(analyze(mol))
        assert parent["kind"] == kind
        assert (parent.get("stem_en"), parent.get("stem_zh")) == kr.parent_names(kind)


def test_select_parent_preserves_existing_stem(monkeypatch: pytest.MonkeyPatch) -> None:
    from namepredict.layer2 import candidates

    parent = {"kind": "phenol", "stem_en": "custom", "stem_zh": "自定义", "mol": object()}
    monkeypatch.setattr(candidates, "_collect_candidates", lambda _: [parent])
    selected = select_parent({"mol": parent["mol"]})
    assert selected["stem_en"] == "custom"
    assert selected["stem_zh"] == "自定义"


def _p(kind: str, **kw) -> dict:
    return {"kind": kind, **kw}


def test_score_parent_tuple_order() -> None:
    """Legacy score consumes the kind-rank fallback for old-style parents."""
    ox = _score_parent({}, _p("oxolane", chain=[0, 1, 2, 3, 4], n_carbons=4))
    eth = _score_parent({}, _p("ether", chain=[0, 1], n_carbons=2))
    assert eth[:2] == (0, 0) and ox[:2] == (0, 0)
    pyr = _score_parent({}, _p("pyridine", chain=list(range(6)), n_carbons=5))
    alk = _score_parent({}, _p("alkane", chain=[0, 1, 2], n_carbons=3))
    assert pyr[3] == 1 and alk[3] == 0  # hetero bit
    assert pyr > alk
    bz = _score_parent({}, _p("acid", chain=list(range(7)), n_carbons=7))
    al = _score_parent({}, _p("alcohol", chain=[0, 1], n_carbons=2))
    assert bz[0] == 14 and al[0] == 5
    assert bz > al


REGRESSION = [
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", REGRESSION)
def test_behavior_regression(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
