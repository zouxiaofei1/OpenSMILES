"""Unit tests for agent_loop.pi_runner (Mock + real runner via subprocess mock)."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from agent_loop.pi_runner import MockPiRunner, PiResult, PiRunner


def test_mock_runner_succeeds_by_default(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("MOCK_PI_FAIL", raising=False)
    prompt = tmp_path / "prompt.md"
    prompt.write_text("hello", encoding="utf-8")
    result = MockPiRunner().run(prompt, tmp_path, timeout_s=5)
    assert isinstance(result, PiResult)
    assert result.success is True
    assert result.exit_code == 0
    assert "mock" in result.log.lower()


def test_mock_runner_fails_when_env_set(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("MOCK_PI_FAIL", "1")
    prompt = tmp_path / "prompt.md"
    prompt.write_text("hello", encoding="utf-8")
    result = MockPiRunner().run(prompt, tmp_path, timeout_s=5)
    assert result.success is False
    assert result.exit_code != 0


def test_pi_runner_success(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    prompt = tmp_path / "prompt.md"
    prompt.write_text("do work", encoding="utf-8")

    def fake_run(*_a, **_k):
        return subprocess.CompletedProcess(
            args=["pi", "-p"],
            returncode=0,
            stdout="ok\n",
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = PiRunner().run(prompt, tmp_path, timeout_s=10)
    assert result.success is True
    assert result.exit_code == 0
    assert "ok" in result.log


def test_pi_runner_nonzero_exit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    prompt = tmp_path / "prompt.md"
    prompt.write_text("fail", encoding="utf-8")

    def fake_run(*_a, **_k):
        return subprocess.CompletedProcess(
            args=["pi", "-p"],
            returncode=2,
            stdout="",
            stderr="boom",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = PiRunner().run(prompt, tmp_path, timeout_s=10)
    assert result.success is False
    assert result.exit_code == 2
    assert "boom" in result.log


def test_pi_runner_timeout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    prompt = tmp_path / "prompt.md"
    prompt.write_text("slow", encoding="utf-8")

    def fake_run(*_a, **_k):
        raise subprocess.TimeoutExpired(cmd=["pi", "-p"], timeout=1)

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = PiRunner().run(prompt, tmp_path, timeout_s=1)
    assert result.success is False
    assert result.exit_code != 0
    assert "timeout" in result.log.lower()


def test_pi_runner_missing_binary(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    prompt = tmp_path / "prompt.md"
    prompt.write_text("x", encoding="utf-8")

    def fake_run(*_a, **_k):
        raise FileNotFoundError("pi")

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = PiRunner().run(prompt, tmp_path, timeout_s=5)
    assert result.success is False
    assert result.exit_code == 127
    assert "pi not found" in result.log.lower()


def test_pi_runner_missing_prompt(tmp_path: Path):
    missing = tmp_path / "nope.md"
    result = PiRunner().run(missing, tmp_path, timeout_s=5)
    assert result.success is False
    assert result.exit_code == 1
    assert "prompt not found" in result.log.lower()


def test_cycle_prompt_has_placeholders():
    root = Path(__file__).resolve().parents[2]
    text = (root / "agent_loop" / "prompts" / "cycle.md").read_text(encoding="utf-8")
    for key in (
        "{{iter}}",
        "{{dual}}",
        "{{cluster}}",
        "{{skill_path}}",
        "{{constraints}}",
        "{{bench_cmd}}",
    ):
        assert key in text, f"missing placeholder {key}"
