# IUPAC: P-44 / architecture
# Layer: L2,L5
"""ParentKind registry: single source for scoring meta + L5 parent stems."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer2 import kind_registry as kr
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.scoring import _score_parent
from namepredict.namer import SMILESNNamer

# Frozen snapshot of pre-migration layer5.benzene_names._ARENE_FG
_LEGACY_ARENE_FG = {
    "phenol": ("phenol", "苯酚"),
    "aniline": ("aniline", "苯胺"),
    "benzoic": ("benzoic acid", "苯甲酸"),
    "benzaldehyde": ("benzaldehyde", "苯甲醛"),
    "acetophenone": ("acetophenone", "苯乙酮"),
    "benzonitrile": ("benzonitrile", "苯甲腈"),
    "benzoyl_chloride": ("benzoyl chloride", "苯甲酰氯"),
    "pyridine": ("pyridine", "吡啶"),
    "furan": ("furan", "呋喃"),
    "thiophene": ("thiophene", "噻吩"),
    "pyrrole": ("1H-pyrrole", "吡咯"),
    "imidazole": ("1H-imidazole", "咪唑"),
    "pyrazole": ("1H-pyrazole", "吡唑"),
    "oxazole": ("1,3-oxazole", "恶唑"),
    "thiazole": ("1,3-thiazole", "噻唑"),
    "pyrimidine": ("pyrimidine", "嘧啶"),
    "pyrazine": ("pyrazine", "吡嗪"),
    "pyridazine": ("pyridazine", "哒嗪"),
    "naphthalene": ("naphthalene", "萘"),
    "anthracene": ("anthracene", "蒽"),
    "indole": ("1H-indole", "吲哚"),
    "indazole": ("1H-indazole", "1H-吲唑"),
    "benzofuran": ("benzofuran", "苯并呋喃"),
    "benzothiophene": ("1-benzothiophene", "苯并[b]噻吩"),
    "benzothiazole": ("1,3-benzothiazole", "1,3-苯并噻唑"),
    "benzoxazole": ("1,3-benzoxazole", "1,3-苯并噁唑"),
    "benzimidazole": ("1H-benzimidazole", "1H-苯并咪唑"),
    "quinoline": ("quinoline", "喹啉"),
    "quinazoline": ("quinazoline", "喹唑啉"),
    "quinoxaline": ("quinoxaline", "喹喔啉"),
    "isoquinoline": ("isoquinoline", "异喹啉"),
    "aziridine": ("aziridine", "氮杂环丙烷"),
    "oxirane": ("oxirane", "环氧乙烷"),
    "oxolane": ("oxolane", "氧杂环戊烷"),
    "oxane": ("oxane", "氧杂环己烷"),
    "pyrrolidine": ("pyrrolidine", "吡咯烷"),
    "piperidine": ("piperidine", "哌啶"),
    "morpholine": ("morpholine", "吗啉"),
    "piperazine": ("piperazine", "哌嗪"),
    "dioxolane": ("1,3-dioxolane", "1,3-二氧戊环"),
    "dioxane": ("1,4-dioxane", "1,4-二氧六环"),
    "thiolane": ("thiolane", "硫杂环戊烷"),
}


def test_oxolane_meta() -> None:
    m = kr.get("oxolane")
    assert m is not None
    assert m.ring == "hetero"
    assert m.n_rings == 1
    assert m.retained is True
    assert m.fg_rank == 0
    assert kr.parent_names("oxolane") == ("oxolane", "氧杂环戊烷")




def test_cycloalkane_polycarboxylic_meta_has_authoritative_base_stem() -> None:
    assert kr.parent_names("cycloalkane_polycarboxylic") == ("cycloalkane", "环烷烃")


def test_anthracene_meta() -> None:
    m = kr.get("anthracene")
    assert m is not None
    assert m.ring == "carbo"
    assert m.n_rings == 3
    assert m.retained is True
    assert kr.parent_names("anthracene") == ("anthracene", "蒽")


@pytest.mark.parametrize(
    "kind",
    ["quinazoline", "quinoxaline"],
)
def test_fused_hetero_d2(kind: str) -> None:
    m = kr.get(kind)
    assert m is not None
    assert m.ring == "hetero"
    assert m.n_rings == 2
    assert m.retained is True
    assert kr.is_hetero_ring(kind) == 1
    assert kr.n_rings_of(kind) == 2


def test_ether_meta() -> None:
    m = kr.get("ether")
    assert m is not None
    assert m.fg_rank == 0
    assert m.ring == "none"
    assert m.n_rings == 0
    assert kr.fg_rank("ether") == 0
    assert kr.has_principal_fg("ether") == 0


def test_acid_meta() -> None:
    m = kr.get("acid")
    assert m is not None
    assert m.fg_rank == 14
    assert kr.fg_rank("acid") == 14
    assert kr.has_principal_fg("acid") == 1


def test_unknown_kind_negative() -> None:
    assert kr.get("no_such") is None
    assert kr.parent_names("no_such") is None
    assert kr.fg_rank("no_such") == 0
    assert kr.has_principal_fg("no_such") == 0
    assert kr.is_hetero_ring("no_such") == 0
    assert kr.is_carbo_ring("no_such") == 0
    assert kr.n_rings_of("no_such") == 0
    assert kr.retained_bonus("no_such") == 0


@pytest.mark.parametrize("kind,names", list(_LEGACY_ARENE_FG.items()))
def test_parent_names_match_legacy_arene_fg(kind: str, names: tuple[str, str]) -> None:
    assert kr.parent_names(kind) == names


_FIXED_PARENT_CASES = [
    ("Oc1ccccc1", "phenol"),
    ("Nc1ccccc1", "aniline"),
    ("O=C(O)c1ccccc1", "benzoic"),
    ("n1ccccc1", "pyridine"),
    ("c1ccc2ccccc2c1", "naphthalene"),
]


@pytest.mark.parametrize("smiles,kind", _FIXED_PARENT_CASES)
def test_select_parent_preloads_registry_stem(smiles: str, kind: str) -> None:
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


def test_all_legacy_kinds_registered() -> None:
    assert set(_LEGACY_ARENE_FG) <= set(kr.all_kinds())


def _p(kind: str, **kw) -> dict:
    return {"kind": kind, **kw}


def test_score_parent_tuple_order() -> None:
    """Legacy score consumes the registry compatibility projection."""
    ox = _score_parent({}, _p("oxolane", chain=[0, 1, 2, 3, 4], n_carbons=4))
    eth = _score_parent({}, _p("ether", chain=[0, 1], n_carbons=2))
    assert eth[:2] == (0, 0) and ox[:2] == (0, 0)
    pyr = _score_parent({}, _p("pyridine", chain=list(range(6)), n_carbons=5))
    alk = _score_parent({}, _p("alkane", chain=[0, 1, 2], n_carbons=3))
    assert pyr[3] == 1 and alk[3] == 0  # hetero bit
    assert pyr > alk
    bz = _score_parent({}, _p("benzoic", chain=list(range(6)), n_carbons=7))
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
