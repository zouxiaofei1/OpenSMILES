# IUPAC: P-44
# Layer: L2,L3,L5
"""Coverage-gated candidate retry across all ranked P-44 candidates."""
from __future__ import annotations

import time

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer, _run_candidates
from namepredict.types import NameResult


def _result(name: str) -> NameResult:
    return NameResult(en=name, zh=f"中-{name}", success=True, source="iupac", time_ms=0.0, meta={})


def _install_phases(monkeypatch: pytest.MonkeyPatch, complete: dict[str, bool]) -> list[str]:
    from namepredict import namer

    candidates = [
        {"kind": kind, "chain": [0], "principal_group_count": 1}
        for kind in complete
    ]
    calls: list[str] = []
    monkeypatch.setattr(namer, "iter_parent_candidates", lambda info: candidates)
    monkeypatch.setattr(
        namer, "_prepare_candidate",
        lambda info, parent, **kw: (parent, [], complete[parent["kind"]]),
    )
    monkeypatch.setattr(
        namer, "_assemble_candidate",
        lambda parent, subst, **kw: calls.append(parent["kind"]) or _result(parent["kind"]),
    )
    return calls


def test_high_phase_partial_keeps_iupac_seniority(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _install_phases(monkeypatch, {"acid": False, "alcohol": True})
    result = _run_candidates({}, depth=0, t0=time.perf_counter())
    assert result.en == "acid"
    assert result.zh == "中-acid"
    assert calls == ["acid"]
    assert result.meta["coverage_complete"] is False


def test_high_phase_complete_wins_without_lower_phase_assembly(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _install_phases(monkeypatch, {"acid": True, "alcohol": True})
    result = _run_candidates({}, depth=0, t0=time.perf_counter())
    assert result.en == "acid"
    assert result.zh == "中-acid"
    assert calls == ["acid"]
    assert result.meta["coverage_complete"] is True


def test_partial_fallback_uses_highest_phase_only_after_complete_scan(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _install_phases(monkeypatch, {"acid": False, "alcohol": False})
    result = _run_candidates({}, depth=0, t0=time.perf_counter())
    assert result.en == "acid"
    assert result.zh == "中-acid"
    assert calls == ["acid"]
    assert result.meta["coverage_complete"] is False
    assert result.meta["fallback"] == "no_coverage_gate"


def test_unassemblable_high_phase_never_capability_downgrades(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = _install_phases(monkeypatch, {"acid": False, "alcohol": True})
    from namepredict import namer

    monkeypatch.setattr(
        namer, "_assemble_candidate",
        lambda parent, subst, **kw: calls.append(parent["kind"]) or None,
    )
    result = _run_candidates({}, depth=0, t0=time.perf_counter())
    assert not result.success
    assert calls == ["acid"]


def test_methylbutylbenzene_prefers_benzene_parent():
    r = SMILESNNamer().name("c1ccc(cc1)CC(C)CC")
    assert r.success
    assert normalize_en(r.en) == normalize_en("(2-methylbutyl)benzene")
    assert normalize_zh(r.zh) == normalize_zh("(2-甲基丁基)苯")
