"""Unit tests for agent_loop.git_gate using a temporary git repository."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from agent_loop.git_gate import GitGate, GitGateError


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


def test_snapshot_returns_sha(tmp_path: Path):
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    sha = gate.snapshot()
    assert len(sha) >= 7
    assert sha == _run_git(repo, "rev-parse", "HEAD")


def test_commit_allowed_path_and_revert(tmp_path: Path):
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    base = gate.snapshot()

    target = repo / "src" / "namepredict" / "layer1" / "rule.py"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("x = 1\n", encoding="utf-8")

    rel = "src/namepredict/layer1/rule.py"
    sha = gate.commit("feat: add rule", [rel])
    assert sha != base
    assert target.read_text(encoding="utf-8") == "x = 1\n"
    assert _run_git(repo, "rev-parse", "HEAD") == sha

    target.write_text("x = 2\n", encoding="utf-8")
    gate.revert(base)
    assert _run_git(repo, "rev-parse", "HEAD") == base
    assert not target.exists()


def test_commit_rejects_disallowed_path(tmp_path: Path):
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    bad = repo / "benchmarks" / "benchmark.py"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text("score = 0\n", encoding="utf-8")
    with pytest.raises(GitGateError):
        gate.commit("bad", ["benchmarks/benchmark.py"])


def test_allowed_paths_includes_safe_prefixes(tmp_path: Path):
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    allowed = gate.allowed_paths
    assert any(p.replace("\\", "/").startswith("src/namepredict") for p in allowed)
    assert any("tests/unit" in p.replace("\\", "/") for p in allowed)
    assert any("agent_loop" in p.replace("\\", "/") for p in allowed)
