"""Git snapshot / commit / allowlist-only revert gate."""

from __future__ import annotations

import subprocess
from pathlib import Path, PurePosixPath

class GitGateError(Exception):
    """Raised when git gate rejects a path or git fails."""


# Safe prefixes the agent loop may commit (posix-style, relative).
DEFAULT_ALLOWED = (
    "src/namepredict/",
    "tests/unit/",
    "agent_loop/memory/",
    "skills/chem-tdd-skill/",
)

# Production paths restored + untracked-cleaned on revert (never wipe memory).
CLEAN_ON_REVERT = (
    "src/namepredict/",
    "tests/unit/",
    "skills/chem-tdd-skill/",
)


class GitGate:
    """Restrict commits to allowed paths; restore only allowlisted paths on revert."""

    def __init__(
        self,
        cwd: Path,
        allowed_paths: tuple[str, ...] | None = None,
    ) -> None:
        self.cwd = Path(cwd)
        self.allowed_paths = tuple(allowed_paths or DEFAULT_ALLOWED)

    def snapshot(self) -> str:
        """Return current HEAD sha."""
        return self._git("rev-parse", "HEAD")

    def commit(self, msg: str, paths: list[str]) -> str:
        """Stage allowed paths and commit; return new HEAD sha."""
        rels = [self._normalize(p) for p in paths]
        self._assert_allowed(rels)
        self._git("add", "--", *rels)
        self._git("commit", "-m", msg, "--", *rels)
        return self.snapshot()

    def revert(self, sha: str) -> None:
        """Restore allowlist paths to sha; clean untracked under CLEAN_ON_REVERT.

        Does not ``git reset --hard`` the whole tree, so dirty tracked files
        outside the allowlist survive.
        """
        for prefix in CLEAN_ON_REVERT:
            self._try_checkout(sha, prefix)
        self._git("clean", "-fd", "--", *CLEAN_ON_REVERT)

    def _try_checkout(self, sha: str, path: str) -> None:
        try:
            self._git("checkout", sha, "--", path)
        except GitGateError:
            pass  # path may not exist at sha (untracked-only; clean handles it)

    def _assert_allowed(self, rels: list[str]) -> None:
        for rel in rels:
            if not self._is_allowed(rel):
                raise GitGateError(f"path not allowed: {rel}")

    def _is_allowed(self, rel: str) -> bool:
        posix = rel.replace("\\", "/")
        for prefix in self.allowed_paths:
            p = prefix.replace("\\", "/")
            if posix == p.rstrip("/") or posix.startswith(p):
                return True
        return False

    def _normalize(self, path: str) -> str:
        raw = path.replace("\\", "/")
        pure = PurePosixPath(raw)
        if pure.is_absolute() or ".." in pure.parts:
            raise GitGateError(f"invalid path: {path}")
        return str(pure)

    def _git(self, *args: str) -> str:
        try:
            r = self._run_git(args)
        except subprocess.CalledProcessError as exc:
            raise GitGateError(self._err_text(exc)) from exc
        return r.stdout.strip()

    def _run_git(self, args: tuple[str, ...]) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", *args],
            cwd=self.cwd,
            check=True,
            capture_output=True,
            text=True,
        )

    @staticmethod
    def _err_text(exc: subprocess.CalledProcessError) -> str:
        return (exc.stderr or exc.stdout or str(exc)).strip()
