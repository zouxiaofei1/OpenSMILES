"""Dev server with auto-restart on src/server changes.

Why not `uvicorn --reload`?  On Windows the uvicorn reloader restarts the
spawned worker by sending CTRL_C_EVENT via `os.kill(pid, CTRL_C_EVENT)`, which
under the hood is GenerateConsoleCtrlEvent(CTRL_C_EVENT, pid).  That call only
works when `pid` is a *process-group* ID, and a multiprocessing-spawned uvicorn
worker is not a group leader - so the call fails (WinError 87), the worker never
exits, the reloader's join() hangs, and the new code never loads.  (Verified on
Windows 11 + Python 3.13 + uvicorn 0.51.)  watchfiles' `watch()` also fails here
(RustNotify: "Input watch path is neither a file nor a directory"), so we poll
mtime snapshots instead and restart uvicorn on change, killing the whole process
tree so no stale child keeps holding the port.

Run from the project root:
    .venv\\Scripts\\python.exe server\\backend\\dev_server.py
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# 8666 rather than the old 8766: 8766 falls inside a Windows excluded port range
# (8726-8825) reserved by Hyper-V/WSL/Docker, so bind() fails with WinError 10013
# no matter what is or is not listening.  Check the current reservations with
# `netsh interface ipv4 show excludedportrange protocol=tcp` before moving this.
PORT = 8666
UVICORN = [
    str(ROOT / ".venv" / "Scripts" / "uvicorn.exe"),
    "server.backend.app:app",
    "--host", "127.0.0.1",
    "--port", str(PORT),
]
SKIP_DIRS = {
    "__pycache__", ".git", ".venv", "node_modules",
    ".mypy_cache", ".pytest_cache", ".ruff_cache",
}


def kill_port_tree(port: int) -> None:
    """Kill every process listening on `port`, tree included.

    Covers plain uvicorn (single worker) and a stale reloader/worker split, so a
    leftover child never keeps the port bound.
    """
    try:
        out = subprocess.run(["netstat", "-ano"], capture_output=True).stdout
    except Exception:
        return
    pids = set()
    for line in out.splitlines():
        line = line.decode(errors="ignore")
        if f":{port}" in line and "LISTENING" in line:
            parts = line.split()
            if parts:
                pids.add(parts[-1])
    for pid in pids:
        subprocess.run(["taskkill", "/F", "/T", "/PID", pid], capture_output=True)


def snapshot() -> dict[str, float]:
    """Map of .py file path -> mtime under src/ and server/."""
    snap: dict[str, float] = {}
    for base in ("src", "server"):
        base_dir = ROOT / base
        if not base_dir.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(base_dir):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            for fn in filenames:
                if not fn.endswith(".py"):
                    continue
                p = os.path.join(dirpath, fn)
                try:
                    snap[p] = os.path.getmtime(p)
                except OSError:
                    pass
    return snap


def start_server() -> subprocess.Popen:
    print(f"[dev] uvicorn on http://127.0.0.1:{PORT}/", flush=True)
    return subprocess.Popen(UVICORN, cwd=str(ROOT))


def restart(proc: subprocess.Popen | None) -> subprocess.Popen:
    if proc is not None and proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass
    kill_port_tree(PORT)
    time.sleep(0.5)
    return start_server()


def main() -> None:
    print("[dev] watching src/ and server/ for changes", flush=True)
    proc = start_server()
    prev = snapshot()
    try:
        while True:
            time.sleep(1.0)
            cur = snapshot()
            if cur == prev:
                continue
            changed = [
                p for p in cur
                if p not in prev or cur[p] != prev[p]
            ] or [p for p in prev if p not in cur]
            print("[dev] change detected - restarting uvicorn", flush=True)
            for p in changed:
                print(f"      {p}", flush=True)
            proc = restart(proc)
            prev = cur
    except KeyboardInterrupt:
        print("[dev] stopping", flush=True)
    finally:
        if proc is not None and proc.poll() is None:
            proc.terminate()
        kill_port_tree(PORT)


if __name__ == "__main__":
    main()
