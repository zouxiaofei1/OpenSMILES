"""Layer benchmark API: sample benchmark rows by tier, generate LLM Layer0 ground
truth via data/run_layer0_llm.py, and score the real pipeline's Layer0 output
against it.

Endpoints (prefix /api/v1):
  POST /layer-benchmark/sample         {ratio, layer}   — stratified sample + start LLM gen
  GET  /layer-benchmark/sample-status                    — poll LLM ground-truth progress
  GET  /layer-benchmark/data                             — current sample/truth state
  POST /layer-benchmark/score          {layer}           — score real L0 vs LLM truth
  GET  /layer-benchmark/score-result                     — last score result

Files:
  data/layer0_sample.json     sampled benchmark rows (input to run_layer0_llm)
  data/layer0_output.json     LLM Layer0 ground truth (written by run_layer0_llm)
  data/layer0_progress.json   LLM run progress (done indices, removed on finish)
  data/layer0_score.json      last scoring result
"""

from __future__ import annotations

import json
import random
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field
from rdkit import Chem

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt

router = APIRouter(prefix="/api/v1", tags=["layer-benchmark"])

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "merged_benchmark.json"
SAMPLE = ROOT / "data" / "layer0_sample.json"
LLM_OUT = ROOT / "data" / "layer0_output.json"
LLM_PROGRESS = ROOT / "data" / "layer0_progress.json"
SCORE = ROOT / "data" / "layer0_score.json"

_LLM_RUNNER = ROOT / "data" / "run_layer0_llm.py"
SPAWN_LOG = ROOT / "data" / "layer0_spawn_error.log"
SAMPLE_META = ROOT / "data" / "layer0_sample_meta.json"  # 持久化上次抽样 ratio/n（重启后仍可读）

# 单物质分析复用 data/run_layer0_llm 的 prompt 与解析（单点维护）
_llm_runner = None


def _import_llm_runner():
    """Lazily import run_layer0_llm (data/ is not a package)."""
    global _llm_runner
    if _llm_runner is not None:
        return _llm_runner
    data_dir = str(ROOT / "data")
    if data_dir not in sys.path:
        sys.path.insert(0, data_dir)
    import run_layer0_llm  # noqa: PLC0415
    _llm_runner = run_layer0_llm
    return _llm_runner

# ── subprocess / mutable state ───────────────────────────────────
_proc: subprocess.Popen | None = None
_lock = threading.Lock()
_sample_ratio = 0.0          # last requested ratio (percent)
_sample_n = 0                # number of sampled rows


def _read_json(path: Path, default: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def _source_rows() -> list[dict]:
    rows = _read_json(SOURCE, [])
    return rows if isinstance(rows, list) else []


def _sample_meta() -> dict[str, Any]:
    meta = _read_json(SAMPLE_META, {})
    return meta if isinstance(meta, dict) else {}


def _sample_rows(rows: list[dict], ratio: float) -> list[dict]:
    """Tier-stratified random sample so easy/hard rows both appear (>=1 per tier).

    ratio is a percentage in (0, 100].
    """
    groups: dict[int, list[int]] = {}
    for i, r in enumerate(rows):
        groups.setdefault(int(r.get("tier", 0) or 0), []).append(i)
    chosen: list[int] = []
    for tier in sorted(groups):
        idxs = groups[tier]
        n = min(len(idxs), max(1, round(len(idxs) * ratio / 100.0)))
        chosen.extend(random.sample(idxs, n))
    chosen.sort()
    return [rows[i] for i in chosen]


# ── LLM ground-truth generation (subprocess) ─────────────────────

def _llm_runner_cmd() -> list[str]:
    return [
        sys.executable, str(_LLM_RUNNER),
        "--input", "data/layer0_sample.json",
        "--output", "data/layer0_output.json",
    ]


def _drain(stream) -> None:
    """Drain a subprocess pipe so the child never blocks on a full buffer."""
    try:
        for _ in stream:
            pass
    except Exception:
        pass


def _start_llm_generation() -> None:
    global _proc
    with _lock:
        if _proc is not None and _proc.poll() is None:
            return
    try:
        _proc = subprocess.Popen(
            _llm_runner_cmd(),
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as exc:
        try:
            SPAWN_LOG.write_text(repr(exc), encoding="utf-8")
        except OSError:
            pass
        _proc = None
        return
    for stream_name in ("stdout", "stderr"):
        stream = getattr(_proc, stream_name)
        if stream is not None:
            threading.Thread(target=_drain, args=(stream,), daemon=True).start()


def _llm_truth_state() -> dict[str, Any]:
    """done/total for the LLM ground-truth run (derived from files).

    total 固定为抽样样本数（progress 文件的 count 只是“已完成数”，
    run_layer0_llm 每次保存时都会刷新，且完成后删除）。
    """
    out = _read_json(LLM_OUT, {})
    results = out.get("results", {}) if isinstance(out, dict) else {}
    done = len(results)
    sample = _read_json(SAMPLE, [])
    total = len(sample) if isinstance(sample, list) else 0
    prog = _read_json(LLM_PROGRESS, None)
    running = (done > 0 and done < total) or (prog is not None and done < total)
    return {
        "done": done,
        "total": total,
        "running": running,
        "ready": total > 0 and done >= total,
    }


# ── real Layer0 (deterministic pipeline) ─────────────────────────

def _real_l0(smiles: str) -> dict[str, Any]:
    """Run the deterministic L0 (preprocess + salt dissociate) on one SMILES."""
    mol = preprocess(smiles)
    if mol is None:
        return {"parseable": False}
    organic, salt_meta = dissociate_salt(mol)
    # n_fragments 与 LLM 侧语义对齐：整分子的独立单元数（点分隔）
    n_frags = len(Chem.GetMolFrags(mol))
    has_c = any(a.GetAtomicNum() == 6 for a in mol.GetAtoms())
    return {
        "parseable": True,
        "classification": "organic" if has_c else "inorganic",
        "salt_meta": salt_meta or None,
        "n_fragments": n_frags,
    }


def _real_salt_type(meta: dict | None) -> str | None:
    if not meta:
        return None
    if "metal" in meta:
        return "alkali_metal"
    if "acid_salt" in meta:
        return "hydrochloride"
    return "other"


def _llm_salt_type(salt: Any) -> str | None:
    if not isinstance(salt, dict):
        return None
    t = salt.get("type")
    return t if isinstance(t, str) else None


def _compare_one(real: dict[str, Any], llm: dict[str, Any]) -> dict[str, Any]:
    """Compare deterministic L0 vs LLM L0 for one row (per-field match or None)."""
    out: dict[str, Any] = {"complete": False}

    lp = llm.get("parseable")
    out["parseable"] = bool(real.get("parseable")) == bool(lp) if lp is not None else None

    lc = llm.get("classification")
    out["classification"] = (real.get("classification") == lc) if lc in ("organic", "inorganic") else None

    r_salt, l_salt = _real_salt_type(real.get("salt_meta")), _llm_salt_type(llm.get("salt"))
    out["both_salt"] = (r_salt is not None and l_salt is not None)
    out["salt_present"] = ((r_salt is not None) == (l_salt is not None))
    # 盐类型：两边都有盐 → 比较 type；两边都无盐 → 视为一致；否则 → 不一致
    if out["both_salt"]:
        out["salt_type"] = (r_salt == l_salt)
    elif r_salt is None and l_salt is None:
        out["salt_type"] = True
    else:
        out["salt_type"] = False

    ln = llm.get("n_fragments")
    out["n_fragments"] = (int(real.get("n_fragments", 0)) == int(ln)) if isinstance(ln, int) else None

    out["complete"] = all(v is not None for v in out.values())
    return out


# ── scoring ──────────────────────────────────────────────────────

def _run_score() -> dict[str, Any]:
    """Score all sampled rows: real L0 vs LLM L0. Runs synchronously (fast)."""
    sample = _read_json(SAMPLE, [])
    truth = _read_json(LLM_OUT, {})
    results = truth.get("results", {}) if isinstance(truth, dict) else {}
    if not sample or not results:
        return {"ok": False, "error": "没有可用的抽样 / LLM 基准数据，请先「刷新 benchmark 数据」"}

    counts = {k: 0 for k in ("parseable", "classification", "salt_present", "salt_type", "n_fragments", "complete")}
    n_salt_both = 0
    diffs: list[dict] = []
    for idx, row in enumerate(sample):
        real = _real_l0(str(row.get("smiles") or ""))
        llm = results.get(str(idx)) or {}
        cmp = _compare_one(real, llm)
        for k in counts:
            if k == "salt_type":  # salt_type 只在 both_salt 分支累计（分母为有盐行）
                continue
            if cmp[k] is True:
                counts[k] += 1
        if cmp.get("both_salt"):
            n_salt_both += 1
            if cmp["salt_type"] is True:
                counts["salt_type"] += 1
        if not (cmp["parseable"] and cmp["classification"] and cmp["salt_present"]
                and cmp["salt_type"] and cmp["n_fragments"]):
            diffs.append({
                "idx": idx, "smiles": row.get("smiles"),
                "real": {"parseable": real["parseable"], "classification": real.get("classification"),
                         "salt": _real_salt_type(real.get("salt_meta")), "n_fragments": real.get("n_fragments")},
                "llm": {"parseable": llm.get("parseable"), "classification": llm.get("classification"),
                        "salt": _llm_salt_type(llm.get("salt")), "n_fragments": llm.get("n_fragments")},
            })

    n = len(sample)
    pct = lambda c: round(100.0 * c / n, 1) if n else 0.0
    score = {
        "ok": True, "n": n,
        "ratio": _sample_meta().get("ratio", _sample_ratio),
        "parseable": {"ok": counts["parseable"], "n": n, "pct": pct(counts["parseable"])},
        "classification": {"ok": counts["classification"], "n": n, "pct": pct(counts["classification"])},
        "salt_present": {"ok": counts["salt_present"], "n": n, "pct": pct(counts["salt_present"])},
        "salt_type": {"ok": counts["salt_type"], "n": n_salt_both, "pct": pct(counts["salt_type"]) if n_salt_both else None},
        "n_fragments": {"ok": counts["n_fragments"], "n": n, "pct": pct(counts["n_fragments"])},
        "complete": {"ok": counts["complete"], "n": n, "pct": pct(counts["complete"])},
        "diffs": diffs[:50],
    }
    SCORE.write_text(json.dumps(score, ensure_ascii=False, indent=1), encoding="utf-8")
    return score


# ── request models ───────────────────────────────────────────────

class SampleRequest(BaseModel):
    ratio: float = Field(default=20.0, ge=1, le=100, description="抽样比例（百分比）")
    layer: int = Field(default=0, ge=0, le=5)


class ScoreRequest(BaseModel):
    layer: int = Field(default=0, ge=0, le=5)


class AnalyzeOneRequest(BaseModel):
    smiles: str = Field(..., min_length=1, description="单物质 SMILES")


# ── endpoints ────────────────────────────────────────────────────

@router.post("/layer-benchmark/sample")
def layer_sample(req: SampleRequest) -> dict[str, Any]:
    """Stratified sample of merged_benchmark.json, then start LLM L0 generation."""
    global _sample_ratio, _sample_n
    if req.layer != 0:
        return {"ok": False, "error": f"Layer {req.layer} 功能未实装"}
    rows = _source_rows()
    if not rows:
        return {"ok": False, "error": f"数据源为空: {SOURCE}"}
    sample = _sample_rows(rows, req.ratio)
    SAMPLE.write_text(json.dumps(sample, ensure_ascii=False), encoding="utf-8")
    SAMPLE_META.write_text(json.dumps({"ratio": req.ratio, "n": len(sample)}), encoding="utf-8")
    _sample_ratio = req.ratio
    _sample_n = len(sample)
    _start_llm_generation()
    return {"ok": True, "started": True, "n": len(sample), "ratio": req.ratio}


@router.get("/layer-benchmark/sample-status")
def layer_sample_status() -> dict[str, Any]:
    """Poll LLM ground-truth generation progress."""
    state = _llm_truth_state()
    return {"ok": True, **state}


@router.get("/layer-benchmark/data")
def layer_data() -> dict[str, Any]:
    """Current sample / ground-truth / last-score state for the page."""
    sample = _read_json(SAMPLE, [])
    truth = _llm_truth_state()
    score = _read_json(SCORE, None)
    return {
        "ok": True,
        "ratio": _sample_meta().get("ratio", _sample_ratio),
        "sample_n": len(sample) if isinstance(sample, list) else 0,
        "truth": truth,
        "score": score if isinstance(score, dict) else None,
    }


@router.post("/layer-benchmark/score")
def layer_score(req: ScoreRequest) -> dict[str, Any]:
    """Score real Layer0 output vs LLM Layer0 ground truth."""
    if req.layer != 0:
        return {"ok": False, "error": f"Layer {req.layer} 功能未实装"}
    return _run_score()


@router.get("/layer-benchmark/score-result")
def layer_score_result() -> dict[str, Any]:
    """Return the last scoring result."""
    score = _read_json(SCORE, None)
    if isinstance(score, dict):
        return {"ok": True, "score": score}
    return {"ok": False, "error": "还没有跑分结果"}


@router.post("/layer-benchmark/analyze-one")
def layer_analyze_one(req: AnalyzeOneRequest) -> dict[str, Any]:
    """Single-molecule Layer0 analysis via LLM（等价 run_layer0_llm --smiles）."""
    try:
        runner = _import_llm_runner()
        raw = runner._call_llm(f"SMILES: {req.smiles}")
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": f"LLM runner 不可用: {exc}"}
    if raw is None:
        return {"ok": False, "error": "LLM 无有效返回（API 失败或重试耗尽）"}
    entry = runner._sanitize(raw, req.smiles)
    for key in ("id", "error"):
        entry.pop(key, None)
    return {"ok": True, "result": entry}
