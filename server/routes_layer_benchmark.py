"""Layer benchmark API: sample benchmark rows by tier, generate LLM per-layer ground
truth via data/run_layer{layer}_llm.py, and score the real pipeline's layer output
against it. Currently supports layer 0 (preprocessor) and 1 (structural analyzer).

Endpoints (prefix /api/v1, layer in {0,1}):
  POST /layer-benchmark/sample         {ratio, layer} — stratified sample + start LLM gen
  GET  /layer-benchmark/sample-status?layer=N          — poll LLM ground-truth progress
  GET  /layer-benchmark/data?layer=N                   — current sample/truth state
  POST /layer-benchmark/score          {layer}         — score real layer vs LLM truth
  GET  /layer-benchmark/score-result?layer=N           — last score result
  POST /layer-benchmark/analyze-one    {smiles, layer} — single-molecule LLM analysis

Per-layer files (layer = 0 or 1):
  data/layer{N}_sample.json     sampled rows (input to run_layer{N}_llm)
  data/layer{N}_output.json     LLM ground truth (written by runner)
  data/layer{N}_progress.json   LLM run progress (done indices, removed on finish)
  data/layer{N}_score.json      last scoring result
  data/layer{N}_sample_meta.json persisted last ratio/n
"""

from __future__ import annotations

import importlib
import json
import random
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any, Callable

from fastapi import APIRouter, Query
from pydantic import BaseModel, Field
from rdkit import Chem

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt
from namepredict.layer1.analyzer import analyze

router = APIRouter(prefix="/api/v1", tags=["layer-benchmark"])

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "merged_benchmark.json"
SUPPORTED_LAYERS = (0, 1)


def _layer_files(layer: int) -> dict[str, Path]:
    return {
        "sample": ROOT / "data" / f"layer{layer}_sample.json",
        "output": ROOT / "data" / f"layer{layer}_output.json",
        "progress": ROOT / "data" / f"layer{layer}_progress.json",
        "meta": ROOT / "data" / f"layer{layer}_sample_meta.json",
        "score": ROOT / "data" / f"layer{layer}_score.json",
        "runner": ROOT / "data" / f"run_layer{layer}_llm.py",
        "spawn_log": ROOT / "data" / f"layer{layer}_spawn_error.log",
    }


# ── LLM runner 导入（按层缓存）───────────────────────────────────
_llm_runners: dict[int, Any] = {}


def _import_llm_runner(layer: int) -> Any:
    if layer in _llm_runners:
        return _llm_runners[layer]
    data_dir = str(ROOT / "data")
    if data_dir not in sys.path:
        sys.path.insert(0, data_dir)
    mod = importlib.import_module(f"run_layer{layer}_llm")
    _llm_runners[layer] = mod
    return mod


# ── subprocess / mutable state ───────────────────────────────────
_procs: dict[int, subprocess.Popen] = {}
_lock = threading.Lock()
_sample_ratio = 0.0
_sample_n = 0


def _read_json(path: Path, default: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def _source_rows() -> list[dict]:
    rows = _read_json(SOURCE, [])
    return rows if isinstance(rows, list) else []


def _sample_meta(layer: int) -> dict[str, Any]:
    meta = _read_json(_layer_files(layer)["meta"], {})
    return meta if isinstance(meta, dict) else {}


def _sample_rows(rows: list[dict], ratio: float) -> list[dict]:
    """Tier-stratified random sample so easy/hard rows both appear (>=1 per tier)."""
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

def _llm_runner_cmd(layer: int) -> list[str]:
    return [
        sys.executable, str(_layer_files(layer)["runner"]),
        "--input", f"data/layer{layer}_sample.json",
        "--output", f"data/layer{layer}_output.json",
    ]


def _drain(stream) -> None:
    try:
        for _ in stream:
            pass
    except Exception:
        pass


def _start_llm_generation(layer: int) -> None:
    with _lock:
        p = _procs.get(layer)
        if p is not None and p.poll() is None:
            return
    try:
        p = subprocess.Popen(
            _llm_runner_cmd(layer),
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        with _lock:
            _procs[layer] = p
    except Exception as exc:
        try:
            _layer_files(layer)["spawn_log"].write_text(repr(exc), encoding="utf-8")
        except OSError:
            pass
        _procs.pop(layer, None)
        return
    for stream_name in ("stdout", "stderr"):
        stream = getattr(p, stream_name)
        if stream is not None:
            threading.Thread(target=_drain, args=(stream,), daemon=True).start()


def _llm_truth_state(layer: int) -> dict[str, Any]:
    """done/total for the LLM ground-truth run (derived from files + subprocess).

    running 的权威信号是子进程是否存活（_procs）——否则刚启动、第一条进度
    尚未落盘时（done=0 且无 progress 文件）会被误判为“未运行”。
    total 固定为抽样样本数（progress 的 count 只是“已完成数”，完成后删除）。
    """
    files = _layer_files(layer)
    out = _read_json(files["output"], {})
    results = out.get("results", {}) if isinstance(out, dict) else {}
    done = len(results)
    sample = _read_json(files["sample"], [])
    total = len(sample) if isinstance(sample, list) else 0
    prog = _read_json(files["progress"], None)

    proc = _procs.get(layer)
    proc_alive = proc is not None and proc.poll() is None
    proc_code = proc.poll() if proc is not None else None

    running = proc_alive or (done > 0 and done < total) or (prog is not None and done < total)
    ready = total > 0 and done >= total

    error = None
    if proc is not None and proc_code is not None and not ready and total > 0:
        error = f"生成进程退出码 {proc_code}" if proc_code != 0 else "生成进程已退出但结果不完整"

    return {
        "done": done,
        "total": total,
        "running": running,
        "ready": ready,
        "error": error,
    }


# ── real Layer0 ──────────────────────────────────────────────────

def _real_l0(smiles: str) -> dict[str, Any]:
    mol = preprocess(smiles)
    if mol is None:
        return {"parseable": False}
    organic, salt_meta = dissociate_salt(mol)
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


_L0_FIELDS = ("parseable", "classification", "salt_present", "salt_type", "n_fragments")


def _compare_l0(real: dict[str, Any], llm: dict[str, Any]) -> dict[str, Any]:
    """Compare deterministic L0 vs LLM L0 for one row (per-field match or None)."""
    out: dict[str, Any] = {"complete": False}

    lp = llm.get("parseable")
    out["parseable"] = bool(real.get("parseable")) == bool(lp) if lp is not None else None

    lc = llm.get("classification")
    out["classification"] = (real.get("classification") == lc) if lc in ("organic", "inorganic") else None

    r_salt, l_salt = _real_salt_type(real.get("salt_meta")), _llm_salt_type(llm.get("salt"))
    out["both_salt"] = (r_salt is not None and l_salt is not None)
    out["salt_present"] = ((r_salt is not None) == (l_salt is not None))
    if out["both_salt"]:
        out["salt_type"] = (r_salt == l_salt)
    elif r_salt is None and l_salt is None:
        out["salt_type"] = True
    else:
        out["salt_type"] = False

    ln = llm.get("n_fragments")
    out["n_fragments"] = (int(real.get("n_fragments", 0)) == int(ln)) if isinstance(ln, int) else None

    out["complete"] = all(out.get(k) is True for k in _L0_FIELDS)
    return out


# ── real Layer1 (structural analyzer) ────────────────────────────

_FG_KEY = {
    "hydroxyls": "hydroxyl", "carboxyls": "carboxyl", "esters": "ester", "amides": "amide",
    "ketones": "ketone", "aldehydes": "aldehyde", "amines": "amine",
    "quaternary_ammoniums": "quaternary_ammonium", "nitriles": "nitrile",
    "double_bonds": "double_bond", "triple_bonds": "triple_bond",
    "acyl_chlorides": "acyl_chloride", "anhydrides": "anhydride", "thiols": "thiol",
    "ethers": "ether", "sulfides": "sulfide", "nitros": "nitro",
    "phosphates": "phosphate", "phosphonics": "phosphonic", "carbamates": "carbamate",
    "carbonates": "carbonate", "sulfoxides": "sulfoxide", "isocyanates": "isocyanate",
    "isothiocyanates": "isothiocyanate", "ureas": "urea", "hydrazines": "hydrazine",
    "guanidines": "guanidine", "sulfonamides": "sulfonamide", "sulfonates": "sulfonate",
    "sulfonyl_chlorides": "sulfonyl_chloride", "sulfonic_acids": "sulfonic_acid",
    "sulfones": "sulfone", "boronics": "boronic",
}


def _real_l1(smiles: str) -> dict[str, Any]:
    mol = preprocess(smiles)
    if mol is None:
        return {"parseable": False}
    info = analyze(mol)
    fgs = sorted({_FG_KEY[k] for k in _FG_KEY if info.get(k)})
    return {
        "parseable": True,
        "functional_groups": fgs,
        "n_rings": int(info.get("n_rings", 0) or 0),
        "n_ring_systems": int(info.get("n_ring_systems", 0) or 0),
        "has_double_bond": bool(info.get("double_bonds")),
        "has_triple_bond": bool(info.get("triple_bonds")),
    }


_L1_FIELDS = ("parseable", "fg", "n_rings", "n_ring_systems", "has_double_bond", "has_triple_bond")


def _compare_l1(real: dict[str, Any], llm: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {"complete": False}

    lp = llm.get("parseable")
    out["parseable"] = bool(real.get("parseable")) == bool(lp) if lp is not None else None

    lfg = llm.get("functional_groups")
    out["fg"] = (sorted(str(x) for x in lfg) == real.get("functional_groups")) if isinstance(lfg, list) else None

    for key in ("n_rings", "n_ring_systems"):
        lv = llm.get(key)
        out[key] = (int(real.get(key, 0)) == int(lv)) if isinstance(lv, (int, float)) and not isinstance(lv, bool) else None

    for key in ("has_double_bond", "has_triple_bond"):
        lv = llm.get(key)
        out[key] = (bool(real.get(key)) == bool(lv)) if isinstance(lv, bool) else None

    out["complete"] = all(out.get(k) is True for k in _L1_FIELDS)
    return out


def _real_layer(layer: int) -> Callable[[str], dict[str, Any]]:
    return {0: _real_l0, 1: _real_l1}[layer]


def _compare_layer(layer: int) -> Callable[[dict, dict], dict]:
    return {0: _compare_l0, 1: _compare_l1}[layer]


# ── scoring ──────────────────────────────────────────────────────

def _run_score(layer: int) -> dict[str, Any]:
    files = _layer_files(layer)
    sample = _read_json(files["sample"], [])
    truth = _read_json(files["output"], {})
    results = truth.get("results", {}) if isinstance(truth, dict) else {}
    if not sample or not results:
        return {"ok": False, "error": "没有可用的抽样 / LLM 基准数据，请先「刷新 benchmark 数据」"}

    real_fn = _real_layer(layer)
    cmp_fn = _compare_layer(layer)
    # 层专属：参与比较的字段 + 展示用的 diff 字段（真实侧/LLM 侧取值键）
    if layer == 0:
        fields = ("parseable", "classification", "salt_present", "salt_type", "n_fragments")
        diff_keys = ("parseable", "classification", "salt", "n_fragments")
        real_side = lambda r: {
            "parseable": r["parseable"], "classification": r.get("classification"),
            "salt": _real_salt_type(r.get("salt_meta")), "n_fragments": r.get("n_fragments"),
        }
        llm_side = lambda l: {
            "parseable": l.get("parseable"), "classification": l.get("classification"),
            "salt": _llm_salt_type(l.get("salt")), "n_fragments": l.get("n_fragments"),
        }
    else:
        fields = ("parseable", "fg", "n_rings", "n_ring_systems", "has_double_bond", "has_triple_bond")
        diff_keys = ("parseable", "functional_groups", "n_rings", "n_ring_systems", "has_double_bond", "has_triple_bond")
        real_side = lambda r: {k: r.get(k) for k in diff_keys}
        llm_side = lambda l: {k: l.get(k) for k in diff_keys}

    counts = {k: 0 for k in (*fields, "complete")}
    n_special = 0  # layer0 有盐行数（salt_type 的分母）
    diffs: list[dict] = []
    for idx, row in enumerate(sample):
        real = real_fn(str(row.get("smiles") or ""))
        llm = results.get(str(idx)) or {}
        cmp = cmp_fn(real, llm)
        for k in counts:
            if k == "salt_type":
                continue
            if cmp[k] is True:
                counts[k] += 1
        if layer == 0:
            if cmp.get("both_salt"):
                n_special += 1
                if cmp["salt_type"] is True:
                    counts["salt_type"] += 1
        if not all(cmp.get(k) for k in fields):
            diffs.append({"idx": idx, "smiles": row.get("smiles"),
                          "real": real_side(real), "llm": llm_side(llm)})

    n = len(sample)
    pct = lambda c: round(100.0 * c / n, 1) if n else 0.0
    score: dict[str, Any] = {"ok": True, "n": n, "ratio": _sample_meta(layer).get("ratio", _sample_ratio)}
    for k in (*fields, "complete"):
        if layer == 0 and k == "salt_type":
            score[k] = {"ok": counts[k], "n": n_special, "pct": pct(counts[k]) if n_special else None}
        else:
            score[k] = {"ok": counts[k], "n": n, "pct": pct(counts[k])}
    score["diffs"] = diffs[:50]
    files["score"].write_text(json.dumps(score, ensure_ascii=False, indent=1), encoding="utf-8")
    return score


# ── request models ───────────────────────────────────────────────

class SampleRequest(BaseModel):
    ratio: float = Field(default=20.0, ge=1, le=100, description="抽样比例（百分比）")
    layer: int = Field(default=0, ge=0, le=5)


class ScoreRequest(BaseModel):
    layer: int = Field(default=0, ge=0, le=5)


class AnalyzeOneRequest(BaseModel):
    smiles: str = Field(..., min_length=1, description="单物质 SMILES")
    layer: int = Field(default=0, ge=0, le=5)


# ── endpoints ────────────────────────────────────────────────────

@router.post("/layer-benchmark/sample")
def layer_sample(req: SampleRequest) -> dict[str, Any]:
    """Stratified sample of merged_benchmark.json, then start LLM generation."""
    global _sample_ratio, _sample_n
    if req.layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {req.layer} 功能未实装"}
    rows = _source_rows()
    if not rows:
        return {"ok": False, "error": f"数据源为空: {SOURCE}"}
    sample = _sample_rows(rows, req.ratio)
    files = _layer_files(req.layer)
    files["sample"].write_text(json.dumps(sample, ensure_ascii=False), encoding="utf-8")
    files["meta"].write_text(json.dumps({"ratio": req.ratio, "n": len(sample)}), encoding="utf-8")
    _sample_ratio = req.ratio
    _sample_n = len(sample)
    _start_llm_generation(req.layer)
    return {"ok": True, "started": True, "n": len(sample), "ratio": req.ratio}


@router.get("/layer-benchmark/sample-status")
def layer_sample_status(layer: int = Query(default=0, ge=0, le=5)) -> dict[str, Any]:
    """Poll LLM ground-truth generation progress."""
    if layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {layer} 功能未实装"}
    state = _llm_truth_state(layer)
    return {"ok": True, **state}


@router.get("/layer-benchmark/data")
def layer_data(layer: int = Query(default=0, ge=0, le=5)) -> dict[str, Any]:
    """Current sample / ground-truth / last-score state for the page."""
    if layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {layer} 功能未实装"}
    files = _layer_files(layer)
    sample = _read_json(files["sample"], [])
    truth = _llm_truth_state(layer)
    score = _read_json(files["score"], None)
    return {
        "ok": True,
        "ratio": _sample_meta(layer).get("ratio", _sample_ratio),
        "sample_n": len(sample) if isinstance(sample, list) else 0,
        "truth": truth,
        "score": score if isinstance(score, dict) else None,
    }


@router.post("/layer-benchmark/score")
def layer_score(req: ScoreRequest) -> dict[str, Any]:
    """Score real layer output vs LLM layer ground truth."""
    if req.layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {req.layer} 功能未实装"}
    return _run_score(req.layer)


@router.get("/layer-benchmark/score-result")
def layer_score_result(layer: int = Query(default=0, ge=0, le=5)) -> dict[str, Any]:
    """Return the last scoring result for a layer."""
    if layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {layer} 功能未实装"}
    score = _read_json(_layer_files(layer)["score"], None)
    if isinstance(score, dict):
        return {"ok": True, "score": score}
    return {"ok": False, "error": "还没有跑分结果"}


@router.post("/layer-benchmark/analyze-one")
def layer_analyze_one(req: AnalyzeOneRequest) -> dict[str, Any]:
    """Single-molecule LLM analysis（等价 run_layer{N}_llm --smiles）."""
    if req.layer not in SUPPORTED_LAYERS:
        return {"ok": False, "error": f"Layer {req.layer} 功能未实装"}
    try:
        runner = _import_llm_runner(req.layer)
        raw = runner._call_llm(f"SMILES: {req.smiles}")
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": f"LLM runner 不可用: {exc}"}
    if raw is None:
        return {"ok": False, "error": "LLM 无有效返回（API 失败或重试耗尽）"}
    entry = runner._sanitize(raw, req.smiles)
    for key in ("id", "error"):
        entry.pop(key, None)
    return {"ok": True, "result": entry}
