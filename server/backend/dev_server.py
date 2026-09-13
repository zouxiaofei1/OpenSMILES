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

import io
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# uvicorn 的启动失败(traceback)与访问日志同时落一份到 tmp/, 事后可查。
LOG_PATH = ROOT / "tmp" / "dev-server.log"
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

_log: io.TextIOWrapper | None = None  # 当前 uvicorn 的日志句柄
_started_at = 0.0  # 本进程起 uvicorn 的时刻, 用于识别"启动即退出"

# 启动后这么快就退出, 基本只有两种原因(src/ 导入失败、端口被占), 提示指个方向。
_FAST_EXIT_SEC = 3.0


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


def _tee(stream, fh) -> None:
    """子进程输出一路写控制台、一路写日志文件。

    重启时日志句柄会被关掉, 而旧进程可能还在吐最后几行, 故写失败直接放弃。
    """
    for line in stream:
        try:
            sys.stdout.write(line)
            sys.stdout.flush()
            fh.write(line)
            fh.flush()
        except (OSError, ValueError):
            return


def start_server() -> subprocess.Popen:
    """起 uvicorn, stdout/stderr 经 tee 同时进控制台与 tmp/dev-server.log。"""
    global _log, _started_at
    print(f"[dev] uvicorn on http://127.0.0.1:{PORT}/  (日志: {LOG_PATH})", flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    _close_log()
    _log = open(LOG_PATH, "a", encoding="utf-8", errors="replace")
    _log.write(f"\n===== {time.strftime('%Y-%m-%d %H:%M:%S')} uvicorn 启动 =====\n")
    _log.flush()
    proc = subprocess.Popen(
        UVICORN,
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
        bufsize=1,  # uvicorn 日志按行出来, 便于 tee 逐行落盘
    )
    threading.Thread(target=_tee, args=(proc.stdout, _log), daemon=True).start()
    _started_at = time.monotonic()
    return proc


def _close_log() -> None:
    """关掉上一轮的日志句柄, 避免重启时句柄泄漏。"""
    global _log
    if _log is not None:
        try:
            _log.close()
        except OSError:
            pass
        _log = None


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
            # uvicorn 自己死了(导入失败 / 端口被占)不会改任何文件, 只靠快照会一直
            # 停在不监听的状态; 这里主动发现并拉起来, 前端就不必等下一次保存。
            if proc.poll() is not None:
                code = proc.returncode
                print(f"[dev] uvicorn 已退出 (code={code}), 重启中", flush=True)
                if time.monotonic() - _started_at < _FAST_EXIT_SEC:
                    print("[dev] 启动即退出: 多为 src/ 导入失败或端口被占, "
                          f"详见上方输出与 {LOG_PATH}", flush=True)
                    time.sleep(1.0)  # 等端口释放 / 避免失败热循环
                proc = restart(proc)
                prev = snapshot()
                continue
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
        _close_log()


if __name__ == "__main__":
    main()
