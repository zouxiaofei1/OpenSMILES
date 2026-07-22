"""Benchmark preview API: run predictions against merged_benchmark.json.

Data source: data/merged_benchmark.json (4070 gold rows).
Cache:       tools/benchmark_pred_preview_data.json.

GET  /api/v1/benchmark-preview        — serve cached rows + generation status
POST /api/v1/benchmark-preview/refresh — start subprocess generation
GET  /api/v1/benchmark-preview/status  — poll generation progress
"""

from __future__ import annotations

import json
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

from fastapi import APIRouter

router = APIRouter(prefix="/api/v1", tags=["benchmark"])

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "merged_benchmark.json"
CACHE = ROOT / "tools" / "benchmark_pred_preview_data.json"

# Subprocess handle for background generation
_proc: subprocess.Popen | None = None
_proc_lock = threading.Lock()
_gen_total = 0
_captured: list[str] = []  # last stdout lines for diagnostics


def _read_cache() -> list[dict]:
    """Read cached rows, or empty list."""
    if not CACHE.is_file():
        return []
    try:
        rows = json.loads(CACHE.read_text(encoding="utf-8"))
        return rows if isinstance(rows, list) else []
    except Exception:
        return []


def _source_total() -> int:
    """Return the number of rows in merged_benchmark.json, or 0."""
    if not SOURCE.is_file():
        return 0
    try:
        rows = json.loads(SOURCE.read_text(encoding="utf-8"))
        return len(rows) if isinstance(rows, list) else 0
    except Exception:
        return 0


@router.get("/benchmark-preview")
def get_benchmark_preview() -> dict[str, Any]:
    """Return benchmark prediction rows from cache + generation status."""
    rows = _read_cache()
    src_total = _source_total()
    generating = _proc is not None and _proc.poll() is None
    gen_done = len(rows)
    gen_total = _gen_total or src_total

    # If the subprocess just finished, do a final read
    if not generating and _proc is not None and _proc.poll() == 0:
        rows = _read_cache()
        gen_done = len(rows)

    return {
        "rows": rows,
        "cached": len(rows),
        "source_total": src_total,
        "stale": len(rows) < src_total,
        "generating": generating,
        "gen_done": gen_done,
        "gen_total": gen_total,
    }


@router.post("/benchmark-preview/refresh")
def refresh_benchmark_preview(force: bool = False) -> dict[str, Any]:
    """Start background subprocess to generate benchmark cache.

    Set force=true to kill any stuck generation and restart.
    """
    global _proc, _gen_total, _captured

    # Check for a running process
    if _proc is not None:
        if _proc.poll() is None:
            if force:
                # Kill the stuck process and clean up
                try:
                    _proc.kill()
                except Exception:
                    pass
                _proc = None
                _captured = []
            else:
                return {
                    "ok": False,
                    "error": "generation already in progress (use force=true to reset)",
                    "total": _gen_total,
                }
        else:
            # Process already dead — clean up the stale reference
            _proc = None
            _captured = []

    src_total = _source_total()
    if src_total == 0:
        return {"ok": False, "error": f"source data not found or empty: {SOURCE}"}

    _gen_total = src_total
    _captured = []

    # Build a small inline script that does the generation
    script = f'''
import json, sys, time
from pathlib import Path
sys.path.insert(0, r"{ROOT / 'src'}")

# Suppress RDKit C++ warnings — they flood stderr and block the pipe buffer
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh
SOURCE = Path(r"{SOURCE}")
CACHE = Path(r"{CACHE}")

def score_pred(pe, pz, ge, gz):
    en_ok = bool(ge) and normalize_en(pe) == normalize_en(ge)
    zh_ok = bool(gz) and normalize_zh(pz) == normalize_zh(gz)
    if ge and gz: dual = en_ok and zh_ok
    elif ge: dual = en_ok
    elif gz: dual = zh_ok
    else: dual = False
    return {{"en_ok": en_ok, "zh_ok": zh_ok, "ok": dual, "ret": bool(pe or pz)}}

rows = json.loads(SOURCE.read_text(encoding="utf-8"))
total = len(rows)
namer = SMILESNNamer()
payload = []
t0 = time.perf_counter()
for i, row in enumerate(rows):
    smi = str(row.get("smiles") or "")
    ge = row.get("english_name") or ""
    gz = row.get("chinese_name") or ""
    try:
        r = namer.name(smi)
        pe, pz = r.en or "", r.zh or ""
    except Exception:
        pe, pz = "", ""
    sc = score_pred(pe, pz, ge, gz)
    payload.append({{"s": smi, "en": pe, "zh": pz, "ge": ge, "gz": gz, **sc}})
    if (i + 1) % 200 == 0 or (i + 1) == total:
        tmp = CACHE.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
        tmp.replace(CACHE)
        elapsed = time.perf_counter() - t0
        print(f"PROGRESS {{i+1}}/{{total}} {{(i+1)/elapsed:.1f}}r/s", flush=True)

elapsed = time.perf_counter() - t0
n_ok = sum(1 for r in payload if r["ok"])
print(f"DONE {{total}} rows {{elapsed:.1f}}s dual_ok={{n_ok}}", flush=True)
'''

    try:
        _proc = subprocess.Popen(
            [sys.executable, "-c", script],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,  # RDKit warnings go to stderr; suppress to avoid pipe blocking
            text=True,
        )
    except Exception as exc:
        _proc = None
        return {"ok": False, "error": str(exc)}

    # Drain stdout in a daemon thread so the pipe buffer never blocks the child
    def _drain() -> None:
        if _proc is None or _proc.stdout is None:
            return
        try:
            for line in _proc.stdout:
                _captured.append(line)
        except Exception:
            pass

    threading.Thread(target=_drain, daemon=True).start()

    return {"ok": True, "started": True, "total": src_total}


@router.get("/benchmark-preview/status")
def benchmark_status() -> dict[str, Any]:
    """Poll generation progress."""
    global _proc
    done = len(_read_cache())
    total = _gen_total or _source_total()

    proc_running = _proc is not None and _proc.poll() is None
    proc_exited = _proc is not None and _proc.poll() is not None

    # If the process exited but cache is incomplete, it failed
    running = proc_running and done < total

    error = None
    if proc_exited:
        if _proc.poll() != 0:
            tail = "".join(_captured[-10:]) if _captured else ""
            error = tail[-500:] if tail else f"exit code {_proc.poll()}"
        # If process exited successfully but done < total, something went wrong
        elif done < total:
            error = "process exited successfully but cache is incomplete"

    # If no process running and cache incomplete, generation is stale/stuck
    stuck = not proc_running and not proc_exited and done < total and done > 0

    return {
        "running": running,
        "done": done,
        "total": total,
        "error": error,
        "stuck": stuck,
    }
