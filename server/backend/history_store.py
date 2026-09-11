"""History-version support: git commits, worktree lifecycle, result cache.

Three data pages (code-analysis / benchmark / call-graph) can be pointed at a
past git commit. Code analysis is computed statically from git objects; the
benchmark and call-graph need the historical code to run, so a worktree is
created lazily, generation runs there, results are cached under
`tools/history_cache/<commit>/`, and the worktree is removed afterwards.

Concurrency: at most one historical generation runs at a time (benchmark or
call-graph, serialised across commits). New requests get a `busy` error and the
frontend auto-retries.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
CACHE_ROOT = ROOT / "tools" / "history_cache"
WT_BASE = ROOT / ".claude" / "worktrees" / "history"
DATA = ROOT / "benchmarks" / "merged_benchmark.json"

# ---------------------------------------------------------------------------
# git commit listing (module-level cache, 5 min TTL)
# ---------------------------------------------------------------------------

_COMMITS_CACHE: tuple[float, list[dict[str, str]]] | None = None
_COMMITS_LOCK = threading.Lock()


def _git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args], cwd=str(ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace",
    )


def _head() -> str | None:
    r = _git("rev-parse", "HEAD")
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else None


def list_commits(force: bool = False) -> dict[str, Any]:
    """HEAD-lineage commits. Returns {ok, head, commits:[{hash,short,date,subject}]}."""
    global _COMMITS_CACHE
    now = time.time()
    with _COMMITS_LOCK:
        if _COMMITS_CACHE and not force and now - _COMMITS_CACHE[0] < 300:
            return {"ok": True, "head": _head(), "commits": _COMMITS_CACHE[1]}
    r = _git("log", "--format=%H%x1f%cI%x1f%s")
    if r.returncode != 0:
        return {"ok": False, "error": (r.stderr or r.stdout)[-500:] or "git log failed"}
    commits: list[dict[str, str]] = []
    for line in r.stdout.splitlines():
        parts = line.split("\x1f", 2)
        if len(parts) != 3:
            continue
        h, date, subject = parts
        commits.append({"hash": h, "short": h[:7], "date": date, "subject": subject})
    with _COMMITS_LOCK:
        _COMMITS_CACHE = (now, commits)
    return {"ok": True, "head": _head(), "commits": commits}


# ---------------------------------------------------------------------------
# commit resolution
# ---------------------------------------------------------------------------

def resolve_commit(commit: str | None) -> str | None:
    """Validate/normalise a ref → full hash. None/empty/HEAD → None (= current).

    Raises ValueError for an unresolvable ref.
    """
    if not commit:
        return None
    c = commit.strip()
    if not c or c in ("HEAD", "current", "当前"):
        return None
    r = _git("rev-parse", "--verify", f"{c}^{{commit}}")
    if r.returncode != 0 or not r.stdout.strip():
        raise ValueError(commit)
    full = r.stdout.strip()
    return None if full == _head() else full


# ---------------------------------------------------------------------------
# result cache (tools/history_cache/<commit>/{kind}.json, atomic writes)
# ---------------------------------------------------------------------------

def cache_dir(commit: str) -> Path:
    return CACHE_ROOT / commit


def cache_path(commit: str, kind: str) -> Path:
    return cache_dir(commit) / f"{kind}.json"


def data_sig() -> str:
    """Signature of the benchmark data file; cache invalidates when it changes."""
    try:
        st = DATA.stat()
        return f"{st.st_size}:{st.st_mtime_ns}"
    except OSError:
        return "missing"


def read_cache(commit: str, kind: str) -> dict | None:
    p = cache_path(commit, kind)
    if not p.is_file():
        return None
    try:
        obj = json.loads(p.read_text(encoding="utf-8"))
        return obj if isinstance(obj, dict) else None
    except Exception:
        return None


def write_cache(commit: str, kind: str, obj: dict) -> None:
    d = cache_dir(commit)
    d.mkdir(parents=True, exist_ok=True)
    p = cache_path(commit, kind)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
    tmp.replace(p)


def write_error(commit: str, kind: str, error: str) -> None:
    """Persist a failed generation's error so status can report it after the job is gone."""
    try:
        d = cache_dir(commit)
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{kind}.err").write_text(error, encoding="utf-8")
    except OSError:
        pass


def read_error(commit: str, kind: str) -> str | None:
    p = cache_dir(commit) / f"{kind}.err"
    if not p.is_file():
        return None
    try:
        return p.read_text(encoding="utf-8")[-500:]
    except Exception:
        return None


# ---------------------------------------------------------------------------
# static git reads (code analysis)
# ---------------------------------------------------------------------------

def list_tree(commit: str, prefix: str) -> list[str]:
    r = _git("ls-tree", "-r", "--name-only", commit, "--", prefix)
    if r.returncode != 0:
        return []
    return [ln.strip() for ln in r.stdout.splitlines() if ln.strip()]


def show_file(commit: str, path: str) -> str | None:
    r = _git("show", f"{commit}:{path}")
    return r.stdout if r.returncode == 0 else None


# ---------------------------------------------------------------------------
# worktree lifecycle (ref-counted)
# ---------------------------------------------------------------------------

_WT_STATE: dict[str, dict[str, Any]] = {}  # path -> {"refs": int, "lock": Lock}


def _wt_path(commit: str) -> Path:
    return WT_BASE / commit[:10]


def _wt_valid(path: Path) -> bool:
    return path.exists() and (path / ".git").exists()


def _rm_tree(path: Path) -> None:
    """Remove a worktree dir: git worktree remove → rmtree → background retry."""
    if not path.exists():
        _git("worktree", "prune")
        return
    _git("worktree", "remove", "--force", str(path))
    if not path.exists():
        _git("worktree", "prune")
        return

    def _retry() -> None:
        for _ in range(5):
            time.sleep(1)
            try:
                shutil.rmtree(path)
                break
            except OSError:
                continue
        _git("worktree", "prune")

    threading.Thread(target=_retry, daemon=True).start()


def ensure_worktree(commit: str) -> Path:
    """Create (if needed) the detached worktree for `commit`; refs += 1."""
    path = _wt_path(commit)
    key = str(path)
    with _WT_STATE.setdefault(key, {"refs": 0, "lock": threading.Lock()})["lock"]:
        if not _wt_valid(path):
            _rm_tree(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            r = _git("worktree", "add", "--detach", str(path), commit)
            if r.returncode != 0:
                raise RuntimeError(f"git worktree add failed: {(r.stderr or r.stdout)[-300:]}")
        _WT_STATE[key]["refs"] += 1
    return path


def release_worktree(commit: str) -> None:
    key = str(_wt_path(commit))
    st = _WT_STATE.get(key)
    if not st:
        return
    with st["lock"]:
        st["refs"] -= 1
        if st["refs"] <= 0:
            st["refs"] = 0
            _rm_tree(_wt_path(commit))


def cleanup_orphans() -> None:
    """Remove leftover history worktrees (startup sweep / manual GC)."""
    active = {str(_wt_path(c)) for c in _WT_STATE if _WT_STATE[c]["refs"] > 0}
    if not WT_BASE.is_dir():
        return
    for d in WT_BASE.iterdir():
        if d.is_dir() and str(d) not in active:
            _rm_tree(d)


# ---------------------------------------------------------------------------
# historical generation job registry (global serial)
# ---------------------------------------------------------------------------

class Job:
    __slots__ = ("commit", "kind", "proc", "total", "captured", "done", "error")

    def __init__(self, commit: str, kind: str) -> None:
        self.commit = commit
        self.kind = kind
        self.proc: subprocess.Popen | None = None
        self.total = 0
        self.captured: list[str] = []
        self.done = 0
        self.error: str | None = None


_HIST_JOBS: dict[tuple[str, str], Job] = {}
_HIST_JOBS_LOCK = threading.Lock()


def _any_running() -> Job | None:
    for j in _HIST_JOBS.values():
        if j.proc is not None and j.proc.poll() is None:
            return j
    return None


def get_job(commit: str, kind: str) -> Job | None:
    with _HIST_JOBS_LOCK:
        return _HIST_JOBS.get((commit, kind))


def register_job(commit: str, kind: str, proc: subprocess.Popen, total: int) -> tuple[bool, str]:
    """Atomically check global busy + register. Returns (ok, error)."""
    with _HIST_JOBS_LOCK:
        busy = _any_running()
        if busy:
            return False, f"busy: {busy.commit[:7]} {busy.kind} 生成中"
        job = Job(commit, kind)
        job.proc = proc
        job.total = total
        _HIST_JOBS[(commit, kind)] = job
    return True, ""


def finish_job(commit: str, kind: str, error: str | None = None) -> Job | None:
    with _HIST_JOBS_LOCK:
        job = _HIST_JOBS.pop((commit, kind), None)
    if job:
        if error:
            write_error(commit, kind, error)
        release_worktree(commit)
    return job
