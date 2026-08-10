"""Benchmark run API: trigger parallel benchmark and return scored results.

POST /api/v1/benchmark-run          — start a benchmark run
GET  /api/v1/benchmark-run/status   — poll progress
GET  /api/v1/benchmark-run/result   — get last completed result
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field

router = APIRouter(prefix="/api/v1", tags=["benchmark-run"])

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_MODULE = "benchmarks.benchmark_parallel"
DATA_DIR = ROOT / "data"
DEFAULT_DATA_FILE = "merged_benchmark.json"


def _data_path(data_file: str | None) -> Path:
    """Resolve a data file relative to the data dir, with traversal guard.

    Accepts a bare filename (resolved against DATA_DIR). Slashes are rejected
    so callers cannot reach arbitrary files elsewhere on disk.
    """
    name = (data_file or DEFAULT_DATA_FILE).strip()
    if "/" in name or "\\" in name or name in ("", ".", ".."):
        raise ValueError("invalid data file name")
    p = (DATA_DIR / name).resolve()
    if p.suffix.lower() != ".json" or not str(p).startswith(str(DATA_DIR.resolve())):
        raise ValueError("invalid data file path")
    return p

# ── mutable state (single-run; new run replaces old) ──────────────
_proc: subprocess.Popen | None = None
_lock = threading.Lock()
_run_total = 0
_run_done = 0
_run_running = False
_run_elapsed = 0.0
_run_rate = 0.0
_run_eta = 0.0
_run_so_far_en_ok = 0
_run_so_far_en_n = 0
_run_so_far_zh_ok = 0
_run_so_far_zh_n = 0
_run_so_far_dual_ok = 0
_run_so_far_dual_n = 0
# Final result from last completed run
_last_result: dict[str, Any] | None = None
_last_error: str | None = None


# ── progress parsing ─────────────────────────────────────────────

# "progress [#########-----] 1500/4070 (36.9%)  86.8 row/s  elapsed=17.3s  eta=29.6s"
_PROGRESS_RE = re.compile(
    r"progress\s+\[.*?\]\s+(\d+)/(\d+)\s+\([\d.]+%\)\s+([\d.]+)\s+row/s\s+elapsed=([\d.]+)s\s+eta=([\d.]+)s"
)
# "  so far  EN=500/1500 (33.3%)  ZH=200/1500 (13.3%)  dual=180/1500 (12.0%)"
_SOFAR_RE = re.compile(
    r"so far\s+EN=(\d+)/(\d+)\s+.*?\s+ZH=(\d+)/(\d+)\s+.*?\s+dual=(\d+)/(\d+)"
)
# "en=36.9% (1500/4070) zh=24.6% (1000/4070) dual=15.5% (629/4070) fails=3441"
_SUMMARY_RE = re.compile(
    r"^en=([\d.]+)%\s+\((\d+)/(\d+)\)\s+zh=([\d.]+)%\s+\((\d+)/(\d+)\)\s+dual=([\d.]+)%\s+\((\d+)/(\d+)\)\s+fails=(\d+)"
)
# "workers=16 elapsed=58.92s avg=14.48ms/row (n=4070)"  (printed on stderr with --time)
_WORKERS_RE = re.compile(
    r"workers=(\d+)\s+elapsed=([\d.]+)s\s+avg=([\d.]+)ms/row\s+\(n=(\d+)\)"
)


def _parse_progress(line: str) -> dict[str, Any] | None:
    m = _PROGRESS_RE.search(line)
    if not m:
        return None
    return {
        "done": int(m.group(1)),
        "total": int(m.group(2)),
        "rate": float(m.group(3)),
        "elapsed": float(m.group(4)),
        "eta": float(m.group(5)),
    }


def _parse_sofar(line: str) -> dict[str, Any] | None:
    m = _SOFAR_RE.search(line)
    if not m:
        return None
    return {
        "en_ok": int(m.group(1)), "en_n": int(m.group(2)),
        "zh_ok": int(m.group(3)), "zh_n": int(m.group(4)),
        "dual_ok": int(m.group(5)), "dual_n": int(m.group(6)),
    }


def _parse_summary(line: str) -> dict[str, Any] | None:
    m = _SUMMARY_RE.match(line)
    if not m:
        return None
    return {
        "acc_en": float(m.group(1)),
        "ok_en": int(m.group(2)), "n_en": int(m.group(3)),
        "acc_zh": float(m.group(4)),
        "ok_zh": int(m.group(5)), "n_zh": int(m.group(6)),
        "acc_dual": float(m.group(7)),
        "ok_dual": int(m.group(8)), "n_dual": int(m.group(9)),
        "fails": int(m.group(10)),
    }


def _parse_workers(line: str) -> dict[str, Any] | None:
    m = _WORKERS_RE.search(line)
    if not m:
        return None
    return {
        "workers": int(m.group(1)),
        "elapsed_sec": float(m.group(2)),
        "avg_ms_per_row": float(m.group(3)),
        "n": int(m.group(4)),
    }


# ── runner ────────────────────────────────────────────────────────

def _run_benchmark(data_path: Path, workers: int, timeout: float,
                   limit: int | None) -> None:
    """Run benchmark_parallel in a subprocess; update globals with progress."""
    global _proc, _run_total, _run_done, _run_running, _run_elapsed, _run_rate, _run_eta
    global _run_so_far_en_ok, _run_so_far_en_n, _run_so_far_zh_ok, _run_so_far_zh_n
    global _run_so_far_dual_ok, _run_so_far_dual_n, _last_result, _last_error

    with _lock:
        _last_error = None
        _last_result = None
        _run_total = 0
        _run_done = 0
        _run_running = True
        _run_elapsed = 0.0
        _run_rate = 0.0
        _run_eta = 0.0

    cmd = [
        sys.executable, "-m", BENCHMARK_MODULE,
        "--data", str(data_path),
        "--timeout", str(timeout),
        "--workers", str(workers),
        "--no-snapshot",
        "--time",
    ]
    if limit is not None:
        cmd += ["--limit", str(limit)]

    try:
        with _lock:
            _proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
    except Exception as exc:
        with _lock:
            _run_running = False
            _last_error = str(exc)
        return

    # Drain stdout + stderr in threads
    stdout_lines: list[str] = []
    stderr_lines: list[str] = []

    def _drain_out() -> None:
        if _proc is None or _proc.stdout is None:
            return
        try:
            for line in _proc.stdout:
                line = line.rstrip("\n")
                stdout_lines.append(line)
                # Parse progress
                prog = _parse_progress(line)
                if prog:
                    with _lock:
                        _run_done = prog["done"]
                        _run_total = prog["total"]
                        _run_rate = prog["rate"]
                        _run_elapsed = prog["elapsed"]
                        _run_eta = prog["eta"]
                    continue
                sofar = _parse_sofar(line)
                if sofar:
                    with _lock:
                        _run_so_far_en_ok = sofar["en_ok"]
                        _run_so_far_en_n = sofar["en_n"]
                        _run_so_far_zh_ok = sofar["zh_ok"]
                        _run_so_far_zh_n = sofar["zh_n"]
                        _run_so_far_dual_ok = sofar["dual_ok"]
                        _run_so_far_dual_n = sofar["dual_n"]
        except Exception:
            pass

    def _drain_err() -> None:
        if _proc is None or _proc.stderr is None:
            return
        try:
            for line in _proc.stderr:
                stderr_lines.append(line.rstrip("\n"))
        except Exception:
            pass

    t_out = threading.Thread(target=_drain_out, daemon=True)
    t_err = threading.Thread(target=_drain_err, daemon=True)
    t_out.start()
    t_err.start()
    t_out.join()
    t_err.join()

    _proc.wait()

    # Parse final result from stdout
    with _lock:
        _run_running = False
        if _proc.returncode != 0:
            tail = "\n".join(stderr_lines[-20:]) if stderr_lines else ""
            _last_error = tail[-500:] if tail else f"exit code {_proc.returncode}"
            return

        # Build result from parsed summary + timing
        result: dict[str, Any] = {}
        for line in stdout_lines:
            summary = _parse_summary(line)
            if summary:
                result.update(summary)
        for line in stderr_lines:
            timing = _parse_workers(line)
            if timing:
                result.update(timing)
        # Fallback: use running so-far values
        if "n_en" not in result:
            result["ok_en"] = _run_so_far_en_ok
            result["n_en"] = _run_so_far_en_n
            result["ok_zh"] = _run_so_far_zh_ok
            result["n_zh"] = _run_so_far_zh_n
            result["ok_dual"] = _run_so_far_dual_ok
            result["n_dual"] = _run_so_far_dual_n
            result["fails"] = max(0, (_run_so_far_dual_n or 0) - (_run_so_far_dual_ok or 0))
        # Compute accuracy percentages
        for field in ("en", "zh", "dual"):
            ok_key = f"ok_{field}"
            n_key = f"n_{field}"
            acc_key = f"acc_{field}"
            if acc_key not in result:
                n = result.get(n_key, 0)
                ok = result.get(ok_key, 0)
                result[acc_key] = round(100.0 * ok / n, 1) if n else 0.0

        _last_result = result


# ── request model ─────────────────────────────────────────────────

class BenchmarkRunRequest(BaseModel):
    data_file: str | None = Field(
        default=None, description=f"Test data file under data/ (default: {DEFAULT_DATA_FILE})"
    )
    workers: int = Field(default=0, ge=0, le=64, description="Process pool size; 0 = auto (cpu_count-1)")
    timeout: float = Field(default=1.0, ge=0.1, le=60.0, description="Per-row timeout in seconds")
    limit: int | None = Field(default=None, ge=1, description="Optional row limit for quick tests")


# ── endpoints ─────────────────────────────────────────────────────

@router.get("/benchmark-run/datasets")
def benchmark_datasets() -> dict[str, Any]:
    """List available .json test files under data/ (bare filenames)."""
    try:
        files = sorted(
            p.name for p in DATA_DIR.iterdir()
            if p.is_file() and p.suffix.lower() == ".json"
        )
    except OSError as exc:
        return {"ok": False, "error": str(exc)}
    return {"ok": True, "datasets": files}


@router.post("/benchmark-run")
def start_benchmark_run(req: BenchmarkRunRequest) -> dict[str, Any]:
    """Start a parallel benchmark run in a background subprocess."""
    global _proc, _run_running

    if _proc is not None and _proc.poll() is None:
        return {
            "ok": False,
            "error": "A benchmark run is already in progress",
        }

    try:
        data_path = _data_path(req.data_file)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}

    if not data_path.is_file():
        return {"ok": False, "error": f"Data not found: {data_path}"}

    workers = req.workers if req.workers > 0 else max(1, (__import__("os").cpu_count() or 2) - 1)
    timeout = req.timeout
    limit = req.limit

    threading.Thread(
        target=_run_benchmark,
        args=(data_path, workers, timeout, limit),
        daemon=True,
    ).start()

    return {
        "ok": True,
        "started": True,
        "data_file": data_path.name,
        "workers": workers,
        "timeout": timeout,
        "limit": limit,
    }


@router.get("/benchmark-run/status")
def benchmark_run_status() -> dict[str, Any]:
    """Poll current run progress."""
    with _lock:
        running = _run_running
        done = _run_done
        total = _run_total
        rate = _run_rate
        elapsed = _run_elapsed
        eta = _run_eta
        sofar = {
            "en_ok": _run_so_far_en_ok, "en_n": _run_so_far_en_n,
            "zh_ok": _run_so_far_zh_ok, "zh_n": _run_so_far_zh_n,
            "dual_ok": _run_so_far_dual_ok, "dual_n": _run_so_far_dual_n,
        }
        error = _last_error
    return {
        "running": running,
        "done": done,
        "total": total,
        "rate": rate,
        "elapsed": elapsed,
        "eta": eta,
        "sofar": sofar,
        "error": error,
    }


@router.get("/benchmark-run/result")
def benchmark_run_result() -> dict[str, Any]:
    """Return the last completed benchmark result, or empty if none."""
    if _last_result is not None:
        return {"ok": True, "result": _last_result}
    if _last_error is not None:
        return {"ok": False, "error": _last_error}
    return {"ok": False, "error": "No benchmark result yet"}
