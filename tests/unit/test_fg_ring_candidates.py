# IUPAC: P-44.1 / P-41
# Layer: L2
"""Multi-candidate FG/ring producers: class-level tries no longer short-circuit.

Layer2 used to return at most one FG and one ring from internal `or` waterfalls,
so `_FG_RANK` never saw acid vs alcohol vs amine as simultaneous candidates.
This suite requires independent producers, rank arbitration, and no regression
on retained parents (benzoic, tert-butylbenzene, ethanol).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.candidates import (
    _collect_candidates,
    _fg_candidates,
    _ring_candidates,
)
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.scoring import _FG_RANK, _pick_best, _score_parent
from namepredict.namer import SMILESNNamer


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _kinds(cands: list[dict]) -> set[str]:
    return {c["kind"] for c in cands if c is not None}


# --- unit: multi-candidate collection ---

def test_aminoalcohol_emits_both_fg_classes() -> None:
    """NCCO: alcohol and amine must both be candidates (not waterfall-or)."""
    kinds = _kinds(_fg_candidates(_info("NCCO")))
    assert "alcohol" in kinds
    assert "amine" in kinds


def test_hydroxyacid_emits_acid_and_alcohol() -> None:
    kinds = _kinds(_fg_candidates(_info("O=C(O)CCO")))
    assert "acid" in kinds
    assert "alcohol" in kinds


def test_collect_includes_fg_ring_alkane() -> None:
    cands = _collect_candidates(_info("c1ccccc1C(=O)O"))
    kinds = _kinds(cands)
    assert "benzoic" in kinds or "acid" in kinds
    assert "benzene" in kinds
    assert "alkane" in kinds
    assert len(cands) >= 3


def test_ring_candidates_simple_benzene() -> None:
    kinds = _kinds(_ring_candidates(_info("c1ccccc1")))
    assert "benzene" in kinds


# --- unit: scoring rank arbitration ---

def test_score_rank_acid_gt_alcohol_gt_amine_gt_alkane() -> None:
    info = _info("C")
    acid = {"kind": "acid", "chain": [0], "n_carbons": 1}
    alcohol = {"kind": "alcohol", "chain": [0], "n_carbons": 1}
    amine = {"kind": "amine", "chain": [0], "n_carbons": 1}
    alkane = {"kind": "alkane", "chain": [0], "n_carbons": 1}
    sa, so, sm, sk = (
        _score_parent(info, acid),
        _score_parent(info, alcohol),
        _score_parent(info, amine),
        _score_parent(info, alkane),
    )
    assert sa > so > sm > sk
    assert _FG_RANK["acid"] > _FG_RANK["alcohol"] > _FG_RANK["amine"]


def test_pick_best_prefers_alcohol_over_amine() -> None:
    info = _info("NCCO")
    cands = _fg_candidates(info) + [{"kind": "alkane", "chain": [0, 1], "n_carbons": 2}]
    best = _pick_best(info, cands)
    assert best is not None
    assert best["kind"] == "alcohol"


def test_pick_best_prefers_acid_over_alcohol() -> None:
    info = _info("O=C(O)CCO")
    cands = _fg_candidates(info)
    best = _pick_best(info, cands)
    assert best is not None
    assert best["kind"] in ("acid", "alkenoic_acid", "diacid")


# --- integration: retained dual behavior ---

REGRESSION = [
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
    ("c1ccc(cc1)C(C)(C)C", "tert-butylbenzene", "叔丁基苯"),
    ("CCO", "ethanol", "乙醇"),
    ("NCCO", "aminoethanol", "氨基乙醇"),  # alcohol parent + amino prefix
]


@pytest.mark.parametrize("smiles,en,zh", REGRESSION)
def test_regression_names(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_select_parent_kind_on_aminoalcohol() -> None:
    p = select_parent(_info("NCCO"))
    assert p["kind"] == "alcohol"


# --- negative: plain alkane must not become a ring ---

def test_plain_alkane_not_ring() -> None:
    p = select_parent(_info("CCC"))
    assert p["kind"] == "alkane"
    r = SMILESNNamer().name("CCC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("propane")
    assert "cyclo" not in normalize_en(r.en)
