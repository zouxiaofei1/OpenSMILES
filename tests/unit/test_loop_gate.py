"""Gate decision tests for AgentLoop (mock pi + mock bench scores)."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

from agent_loop.config import LoopConfig
from agent_loop.events import EventBus
from agent_loop.git_gate import GitGate
from agent_loop.loop import AgentLoop, CycleResult
from agent_loop.pi_runner import MockPiRunner
from agent_loop.state import StateStore


def _run_git(cwd: Path, *args: str) -> str:
    r = subprocess.run(
        ["git", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return r.stdout.strip()


def _init_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    repo.mkdir()
    _run_git(repo, "init")
    _run_git(repo, "config", "user.email", "test@example.com")
    _run_git(repo, "config", "user.name", "Test User")
    (repo / "README.md").write_text("init\n", encoding="utf-8")
    _run_git(repo, "add", "README.md")
    _run_git(repo, "commit", "-m", "initial")
    return repo


def _report(dual: float, n_fails: int = 10) -> dict[str, Any]:
    fails = [
        {
            "id": f"f{i}",
            "smiles": "C" * (i + 1),
            "features": ["alkane"],
            "en_ok": False,
            "zh_ok": False,
            "source": "tiers",
        }
        for i in range(n_fails)
    ]
    return {
        "acc_dual": dual,
        "acc_en": dual,
        "acc_zh": dual,
        "n_dual": 100,
        "ok_dual": int(dual * 100),
        "fails": fails,
    }


def _make_loop(
    repo: Path,
    scores: list[dict[str, Any]],
    *,
    drop_tol: float = 0.005,
    max_iters: int = 1,
) -> AgentLoop:
    mem = repo / "agent_loop" / "memory"
    mem.mkdir(parents=True, exist_ok=True)
    (repo / "agent_loop" / "prompts").mkdir(parents=True, exist_ok=True)
    tmpl = Path("agent_loop/prompts/cycle.md")
    if tmpl.is_file():
        (repo / "agent_loop" / "prompts" / "cycle.md").write_text(
            tmpl.read_text(encoding="utf-8"), encoding="utf-8"
        )
    else:
        (repo / "agent_loop" / "prompts" / "cycle.md").write_text(
            "{{iter}} {{dual}} {{cluster}} {{skill_path}} {{constraints}} {{bench_cmd}}\n",
            encoding="utf-8",
        )
    cfg = LoopConfig(
        data_path=repo / "data" / "merged.json",
        cwd=repo,
        max_iters=max_iters,
        drop_tol=drop_tol,
        mock_pi=True,
        bench_limit=5,
        k=5,
        target_dual=0.99,
    )
    scores_iter = iter(scores)

    def bench_fn() -> dict[str, Any]:
        return next(scores_iter)

    return AgentLoop(
        cfg,
        namer_factory=lambda: object(),
        pi_runner=MockPiRunner(),
        bus=EventBus(),
        bench_fn=bench_fn,
        pytest_fn=lambda: True,
        lint_fn=lambda: [],
        git=GitGate(repo),
        store=StateStore(mem),
    )


def test_gate_reverts_on_dual_drop(tmp_path: Path):
    repo = _init_repo(tmp_path)
    loop = _make_loop(repo, [_report(0.5, 10), _report(0.4, 12)])
    result = loop.run_once()
    assert isinstance(result, CycleResult)
    assert result.dual_before == 0.5
    assert result.dual_after == 0.4
    assert result.decision == "revert"
    assert result.session_path.is_file()
    session = json.loads(result.session_path.read_text(encoding="utf-8"))
    assert session["decision"] == "revert"


def test_gate_commits_on_dual_up(tmp_path: Path):
    repo = _init_repo(tmp_path)
    # Simulate an allowlisted change so GitGate has something to commit.
    target = repo / "src" / "namepredict" / "layer5" / "note.py"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("x = 1\n", encoding="utf-8")

    loop = _make_loop(repo, [_report(0.5, 10), _report(0.6, 8)])
    result = loop.run_once()
    assert result.dual_before == 0.5
    assert result.dual_after == 0.6
    assert result.decision == "commit"
    assert result.session_path.is_file()
    session = json.loads(result.session_path.read_text(encoding="utf-8"))
    assert session["decision"] == "commit"
    assert session["status"] == "gate_pass"
