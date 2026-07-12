"""Run pi CLI non-interactively for one agent cycle."""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PiResult:
    """Outcome of one pi invocation."""

    success: bool
    log: str
    exit_code: int


class PiRunner:
    """Invoke `pi -p` with prompt file content as the initial message."""

    def run(self, prompt_path: Path, cwd: Path, timeout_s: int) -> PiResult:
        """Run pi in print mode; never raise on process failure."""
        try:
            return self._invoke(prompt_path, cwd, timeout_s)
        except FileNotFoundError as exc:
            return self._missing_result(exc)
        except subprocess.TimeoutExpired as exc:
            return self._timeout_result(exc)
        except OSError as exc:
            return PiResult(False, f"os error: {exc}", 1)

    def _invoke(self, prompt_path: Path, cwd: Path, timeout_s: int) -> PiResult:
        text = self._read_prompt(prompt_path)
        proc = self._run_pi(text, Path(cwd), timeout_s)
        return self._from_proc(proc)

    @staticmethod
    def _read_prompt(prompt_path: Path) -> str:
        path = Path(prompt_path)
        if not path.is_file():
            raise FileNotFoundError(f"prompt not found: {path}")
        return path.read_text(encoding="utf-8")

    @staticmethod
    def _missing_result(exc: FileNotFoundError) -> PiResult:
        msg = str(exc)
        if "prompt not found" in msg:
            return PiResult(False, msg, 1)
        return PiResult(False, f"pi not found: {exc}", 127)

    def _run_pi(
        self, prompt: str, cwd: Path, timeout_s: int
    ) -> subprocess.CompletedProcess[str]:
        # pi -p: print response and exit (non-interactive).
        return subprocess.run(
            ["pi", "-p", prompt],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )

    def _from_proc(self, proc: subprocess.CompletedProcess[str]) -> PiResult:
        log = self._combine(proc.stdout, proc.stderr)
        ok = proc.returncode == 0
        return PiResult(ok, log, int(proc.returncode))

    @staticmethod
    def _combine(stdout: str | None, stderr: str | None) -> str:
        parts = [p for p in (stdout or "", stderr or "") if p]
        return "\n".join(parts).strip()

    @staticmethod
    def _timeout_result(exc: subprocess.TimeoutExpired) -> PiResult:
        partial = ""
        if exc.stdout:
            partial = exc.stdout if isinstance(exc.stdout, str) else exc.stdout.decode()
        msg = f"timeout after {exc.timeout}s"
        log = f"{partial}\n{msg}".strip() if partial else msg
        return PiResult(False, log, 124)


class MockPiRunner:
    """No-op runner for tests; fail when MOCK_PI_FAIL is truthy."""

    def run(self, prompt_path: Path, cwd: Path, timeout_s: int) -> PiResult:
        """Return success unless MOCK_PI_FAIL is set."""
        _ = (prompt_path, cwd, timeout_s)
        if self._should_fail():
            return PiResult(False, "mock pi failed (MOCK_PI_FAIL)", 1)
        return PiResult(True, "mock pi ok", 0)

    @staticmethod
    def _should_fail() -> bool:
        val = os.environ.get("MOCK_PI_FAIL", "").strip().lower()
        return val in {"1", "true", "yes", "on"}
