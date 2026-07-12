"""Tests: fresh-process bench path and allowlist-only git gate integration."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any
from unittest.mock import patch

from agent_loop.config import LoopConfig
from agent_loop.loop import AgentLoop, _bench_argv, _bench_subprocess


def test_bench_argv_includes_json_and_limit():
    cfg = LoopConfig(data_path=Path("data/x.json"), bench_limit=7)
    argv = _bench_argv(cfg)
    assert "-m" in argv and "benchmarks.benchmark" in argv
    assert "--json" in argv
    assert "--limit" in argv and "7" in argv
    joined = " ".join(str(a).replace("\\", "/") for a in argv)
    assert "data/x.json" in joined


def test_default_bench_uses_subprocess(tmp_path: Path):
    """Without bench_fn, AgentLoop._bench must call the cold subprocess helper."""
    cfg = LoopConfig(
        data_path=tmp_path / "data.json",
        cwd=tmp_path,
        mock_pi=True,
        bench_limit=1,
        max_iters=1,
    )
    (tmp_path / "agent_loop" / "memory").mkdir(parents=True)
    (tmp_path / "agent_loop" / "prompts").mkdir(parents=True)
    (tmp_path / "agent_loop" / "prompts" / "cycle.md").write_text(
        "{{iter}} {{dual}} {{cluster}} {{skill_path}} {{constraints}} {{bench_cmd}}\n",
        encoding="utf-8",
    )
    fake = {
        "acc_dual": 0.5,
        "acc_en": 0.5,
        "acc_zh": 0.5,
        "n_dual": 2,
        "ok_dual": 1,
        "fails": [],
    }
    loop = AgentLoop(cfg, bench_fn=None, pytest_fn=lambda: True, lint_fn=lambda: [])
    with patch("agent_loop.loop._bench_subprocess", return_value=fake) as mock_sub:
        out = loop._bench()
    assert out == fake
    mock_sub.assert_called_once()


def test_bench_subprocess_parses_json(tmp_path: Path, monkeypatch):
    """_bench_subprocess loads JSON stdout from a cold process."""
    report = {"acc_dual": 0.42, "fails": [], "acc_en": 0.4, "acc_zh": 0.4}

    class _R:
        returncode = 0
        stdout = json.dumps(report)
        stderr = ""

    def _fake_run(*_a, **_k):
        return _R()

    monkeypatch.setattr(subprocess, "run", _fake_run)
    cfg = LoopConfig(data_path=tmp_path / "d.json", bench_limit=3)
    assert _bench_subprocess(cfg, tmp_path)["acc_dual"] == 0.42
