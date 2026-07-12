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


def test_commit_allowed_path(tmp_path: Path):
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


def test_revert_restores_allowlist_keeps_other_dirty(tmp_path: Path):
    """Allowlist dirty is restored; non-allowlist dirty tracked files survive."""
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    base = gate.snapshot()

    # Tracked allowlist file at base
    tracked = repo / "src" / "namepredict" / "tracked.py"
    tracked.parent.mkdir(parents=True, exist_ok=True)
    tracked.write_text("v1\n", encoding="utf-8")
    _run_git(repo, "add", "src/namepredict/tracked.py")
    _run_git(repo, "commit", "-m", "add tracked")
    base = gate.snapshot()

    tracked.write_text("v2 dirty\n", encoding="utf-8")
    (repo / "README.md").write_text("readme dirty\n", encoding="utf-8")
    untracked = repo / "src" / "namepredict" / "new_rule.py"
    untracked.write_text("new\n", encoding="utf-8")

    gate.revert(base)

    assert _run_git(repo, "rev-parse", "HEAD") == base
    assert tracked.read_text(encoding="utf-8") == "v1\n"
    assert not untracked.exists()
    # Non-allowlist dirty must survive (no hard reset)
    assert (repo / "README.md").read_text(encoding="utf-8") == "readme dirty\n"


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


def test_revert_cleans_untracked_allowlisted_keeps_memory(tmp_path: Path):
    repo = _init_repo(tmp_path)
    gate = GitGate(repo)
    base = gate.snapshot()

    prod = repo / "src" / "namepredict" / "new_rule.py"
    prod.parent.mkdir(parents=True, exist_ok=True)
    prod.write_text("x = 1\n", encoding="utf-8")
    unit = repo / "tests" / "unit" / "test_new.py"
    unit.parent.mkdir(parents=True, exist_ok=True)
    unit.write_text("def test_x():\n    assert True\n", encoding="utf-8")
    mem = repo / "agent_loop" / "memory" / "sessions" / "keep.json"
    mem.parent.mkdir(parents=True, exist_ok=True)
    mem.write_text("{}\n", encoding="utf-8")

    gate.revert(base)
    assert not prod.exists()
    assert not unit.exists()
    assert mem.is_file()
