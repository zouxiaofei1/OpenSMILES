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
    real_try = namer_mod.try_candidate

    def tracking_try(info, parent, *, depth=0, t0=None):
        seen.append(parent.get("kind") or "?")
        return real_try(info, parent, depth=depth, t0=t0)

    monkeypatch.setattr(namer_mod, "try_candidate", tracking_try)
    r = SMILESNNamer().name("c1ccc(cc1)CC(C)CC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("(2-methylbutyl)benzene")
    assert "alkane" in seen
    assert "benzene" in seen
    assert seen.index("alkane") < seen.index("benzene")


def test_all_candidates_fail_returns_failure(monkeypatch):
    mol = Chem.MolFromSmiles("CCO")
    bad1 = {"chain": [0], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({0})}
    bad2 = {"chain": [1], "n_carbons": 1, "kind": "alkane", "owned_atoms": frozenset({1})}

    import namepredict.namer as namer_mod

    monkeypatch.setattr(namer_mod, "iter_parent_candidates", lambda _info: [bad1, bad2])
    # Force incomplete coverage: no named side atoms
    monkeypatch.setattr(namer_mod, "extract_substituents", lambda _info, _p: [])
    r = namer_mod._name_mol(mol, depth=0)
    assert r.success is False
    assert r.en == ""


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
