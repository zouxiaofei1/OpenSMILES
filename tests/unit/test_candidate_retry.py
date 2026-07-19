"""Coverage-gated candidate retry: first incomplete candidate yields to next."""
from __future__ import annotations

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer, try_candidate
from namepredict.layer1.analyzer import analyze
from rdkit import Chem


def test_retry_skips_incomplete_first_candidate(monkeypatch):
    """First candidate fails coverage; second succeeds (real molecule)."""
    # (2-methylbutyl)benzene: alkane first is incomplete; benzene second succeeds
    import namepredict.namer as namer_mod

    seen: list[str] = []
    real_prepare = namer_mod._prepare_candidate

    def tracking_prepare(info, parent):
        seen.append(parent.get("kind") or "?")
        return real_prepare(info, parent)

    monkeypatch.setattr(namer_mod, "_prepare_candidate", tracking_prepare)
    r = SMILESNNamer().name("c1ccc(cc1)CC(C)CC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("(2-methylbutyl)benzene")
    assert "alkane" in seen
    assert "benzene" in seen
    assert seen.index("alkane") < seen.index("benzene")


def test_incomplete_candidates_fall_back_to_partial_name(monkeypatch):
    mol = Chem.MolFromSmiles("CCO")
    bad1 = {"chain": [0], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({0})}
    bad2 = {"chain": [1], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({1})}

    import namepredict.namer as namer_mod

    monkeypatch.setattr(namer_mod, "iter_parent_candidates", lambda _info: [bad1, bad2])
    # Force incomplete coverage: no named side atoms
    monkeypatch.setattr(namer_mod, "extract_substituents", lambda _info, _p: [])
    r = namer_mod._name_mol(mol, depth=0)
    assert r.success is True
    assert r.en == "methane"
    assert r.meta["coverage_complete"] is False


def test_fallback_reuses_first_pass_extraction(monkeypatch):
    """Fallback assembly must not extract the same candidates a second time."""
    mol = Chem.MolFromSmiles("CCO")
    candidates = [
        {"chain": [0], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({0})},
        {"chain": [1], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({1})},
    ]

    import namepredict.namer as namer_mod

    calls: list[int] = []

    def track_extraction(_info, parent):
        calls.append(parent["chain"][0])
        return []

    monkeypatch.setattr(namer_mod, "iter_parent_candidates", lambda _info: candidates)
    monkeypatch.setattr(namer_mod, "extract_substituents", track_extraction)

    namer_mod._name_mol(mol, depth=0)

    assert calls == [0, 1]


def test_methylbutylbenzene_prefers_benzene_parent():
    r = SMILESNNamer().name("c1ccc(cc1)CC(C)CC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("(2-methylbutyl)benzene")
    assert normalize_zh(r.zh) == normalize_zh("(2-甲基丁基)苯")


def test_try_candidate_rejects_gap(monkeypatch):
    mol = Chem.MolFromSmiles("CCO")
    info = analyze(mol)
    parent = {
        "chain": [0],
        "n_carbons": 1,
        "kind": "alkane",
        "owned_atoms": frozenset({0}),
        "mol": mol,
    }
    import namepredict.namer as namer_mod

    monkeypatch.setattr(namer_mod, "extract_substituents", lambda _info, _p: [])
    r = try_candidate(info, parent, depth=0)
    assert r is None
