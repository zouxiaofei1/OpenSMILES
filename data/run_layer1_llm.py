"""LLM-driven Layer1 (Structural Analyzer) 分析器：对样本调用 OpenAI-compatible API，
让 LLM 产出 layer1 结果（官能团 / 环 / 不饱和度），用于与真实管线 layer1 对比。

结构与 run_layer0_llm.py 相同：串行处理、断点续跑、每 SAVE_EVERY 条写盘一次。

用法（在项目根目录）:
    ./.venv/Scripts/python.exe data/run_layer1_llm.py --input data/layer1_sample.json --output data/layer1_output.json
    ./.venv/Scripts/python.exe data/run_layer1_llm.py --smiles 'CCO'   # 单物质，只打印

环境变量：OPENAI_API_KEY / OPENAI_BASE_URL / OPENAI_MODEL / LLM_DELAY（同 layer0）
"""

from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import httpx

# 强制 UTF-8 输出，避免 Windows 控制台把中文按 GBK 解码成乱码
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")

# =============================================================================
# CONFIGURATION
# =============================================================================
_CFG_DEFAULTS = {
    "model": "deepseek-v4-pro",
    "base_url": "https://api.deepseek.com",
    "api_key": "",
    "concurrency": 1,
    "delay": 0.3,
}
SAVE_EVERY      = 5
TIMEOUT         = 60.0
MAX_TOKENS      = 4000
CALL_RETRIES    = 3

ROOT        = Path(__file__).resolve().parents[1]
SETTINGS_PATH = ROOT / "server" / "settings.json"
INPUT_PATH  = ROOT / "data" / "layer1_sample.json"
OUTPUT_PATH = ROOT / "data" / "layer1_output.json"
PROGRESS_PATH = ROOT / "data" / "layer1_progress.json"


def _load_llm_config() -> dict[str, Any]:
    """每次调用读一次 server/settings.json（联动设置页）；env 优先级最高。"""
    cfg = dict(_CFG_DEFAULTS)
    if SETTINGS_PATH.is_file():
        try:
            data = json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            data = {}
        if isinstance(data, dict):
            if data.get("model"):
                cfg["model"] = data["model"]
            if data.get("base_url"):
                cfg["base_url"] = data["base_url"]
            if data.get("api_key"):
                cfg["api_key"] = data["api_key"]
            try:
                cfg["concurrency"] = max(1, int(data.get("concurrency") or 1))
            except (TypeError, ValueError):
                pass
    cfg["model"] = os.environ.get("OPENAI_MODEL", cfg["model"])
    cfg["base_url"] = os.environ.get("OPENAI_BASE_URL", cfg["base_url"])
    cfg["api_key"] = os.environ.get("OPENAI_API_KEY", cfg["api_key"])
    try:
        cfg["concurrency"] = max(1, int(os.environ.get("LLM_CONCURRENCY", cfg["concurrency"])))
    except ValueError:
        pass
    cfg["delay"] = float(os.environ.get("LLM_DELAY", cfg["delay"]))
    return cfg

# 允许的官能团类型（与 server 真实侧 _FG_KEY 派生保持一致）
_FG_ALLOWED = [
    "aldehyde", "amine", "quaternary_ammonium", "nitrile", "double_bond", "triple_bond",
    "acyl_chloride", "anhydride", "thiol", "ether", "sulfide", "nitro",
    "phosphate", "phosphonic", "carbamate", "carbonate", "sulfoxide",
    "isocyanate", "isothiocyanate", "urea", "hydrazine", "guanidine",
    "sulfonamide", "sulfonate", "sulfonyl_chloride", "sulfonic_acid", "sulfone", "boronic",
    "hydroxyl", "carboxyl", "ester", "amide", "ketone",
]
_FG_ALLOWED_STR = ", ".join(_FG_ALLOWED)

SYSTEM_PROMPT = f"""\
You are Layer 1 (Structural Analyzer) of a SMILES→IUPAC naming pipeline. Given a SMILES \
and optional gold IUPAC names, analyze the molecular structure and output:

1. parseable — can the SMILES be parsed into a valid molecule? (true/false)
2. functional_groups — array of functional group types detected. Every name MUST be \
   chosen from this exact list: {_FG_ALLOWED_STR}
   (only include types actually present; if none, return an empty array [])
3. n_rings — total number of rings (0 for acyclic)
4. n_ring_systems — number of independent ring systems (fused/bridged/spiro rings count \
   as one system)
5. has_double_bond — true if any non-aromatic C=C double bond
6. has_triple_bond — true if any C≡C triple bond (not C≡N nitrile)
7. notes — one short sentence.

Respond ONLY with a single JSON object, no markdown fences. Example:
{{"parseable": true, "functional_groups": ["hydroxyl", "carboxyl"], "n_rings": 1, \
"n_ring_systems": 1, "has_double_bond": false, "has_triple_bond": false, \
"notes": "Benzoic acid: one aromatic ring, hydroxyl + carboxyl."}}\
"""


def _build_user_message(row: dict[str, Any]) -> str:
    parts = [f"SMILES: {row.get('smiles', '')}"]
    if row.get("english_name"):
        parts.append(f"English IUPAC name (gold): {row['english_name']}")
    if row.get("chinese_name"):
        parts.append(f"Chinese IUPAC name (gold): {row['chinese_name']}")
    return "\n".join(parts)


# =============================================================================
# API 客户端
# =============================================================================
def _chat_completions(messages: list[dict], *, model: str, api_key: str, base_url: str) -> str:
    payload = {"model": model, "messages": messages, "temperature": 0.0, "max_tokens": MAX_TOKENS}
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    url = base_url.rstrip("/") + "/chat/completions"
    for with_format in (True, False):
        p = dict(payload)
        if with_format:
            p["response_format"] = {"type": "json_object"}
        with httpx.Client(timeout=TIMEOUT) as client:
            resp = client.post(url, headers=headers, json=p)
            if resp.status_code == 200:
                return resp.json()["choices"][0]["message"]["content"]
            if not with_format or "response_format" not in resp.text.lower():
                resp.raise_for_status()
    raise RuntimeError("unreachable")


def _call_llm(user_message: str, cfg: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """调用 LLM 并解析 JSON；空响应/解析失败重试 CALL_RETRIES 次。

    cfg 为 None 时每次调用现读 server/settings.json（保证设置页改动即时生效）。
    """
    cfg = cfg or _load_llm_config()
    delay = float(cfg.get("delay", 0.3))
    messages = [{"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": user_message}]
    for attempt in range(CALL_RETRIES):
        content = ""
        try:
            content = _chat_completions(
                messages,
                model=cfg["model"],
                api_key=cfg.get("api_key", ""),
                base_url=cfg.get("base_url", ""),
            ).strip()
            if not content:
                print(f"  [WARN] empty response (attempt {attempt + 1}/{CALL_RETRIES})", file=sys.stderr)
                time.sleep(delay)
                continue
            return json.loads(content)
        except (json.JSONDecodeError, KeyError, IndexError):
            print(f"  [WARN] bad JSON response (attempt {attempt + 1}/{CALL_RETRIES}): {content[:200]!r}", file=sys.stderr)
        except Exception as exc:  # noqa: BLE001
            print(f"  [ERR] API call failed (attempt {attempt + 1}/{CALL_RETRIES}): {exc}", file=sys.stderr)
        time.sleep(delay)
    return None


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _sanitize(raw: dict[str, Any] | None, smiles: str) -> dict[str, Any]:
    """规整 LLM 输出：functional_groups 仅保留允许列表、类型正确。"""
    raw = raw or {}
    fg_raw = raw.get("functional_groups")
    allowed = set(_FG_ALLOWED)
    fg = sorted({str(x).strip() for x in fg_raw if str(x).strip() in allowed}) if isinstance(fg_raw, list) else []
    return {
        "id": None,
        "smiles": smiles,
        "parseable": _as_bool(raw.get("parseable"), True),
        "functional_groups": fg,
        "n_rings": _as_int(raw.get("n_rings")),
        "n_ring_systems": _as_int(raw.get("n_ring_systems")),
        "has_double_bond": _as_bool(raw.get("has_double_bond")),
        "has_triple_bond": _as_bool(raw.get("has_triple_bond")),
        "notes": str(raw.get("notes") or ""),
        "error": not isinstance(raw, dict),
    }


# =============================================================================
# 结果与进度存取
# =============================================================================
def _load_results() -> dict[str, Any]:
    if not OUTPUT_PATH.exists():
        return {"meta": {}, "results": {}}
    try:
        with OUTPUT_PATH.open(encoding="utf-8") as f:
            data = json.load(f)
        return {"meta": data.get("meta", {}), "results": data.get("results", {})}
    except (json.JSONDecodeError, OSError):
        return {"meta": {}, "results": {}}


def _save_results(payload: dict[str, Any]) -> None:
    tmp = OUTPUT_PATH.with_suffix(".json.tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=1)
    tmp.replace(OUTPUT_PATH)


def _load_done() -> set[int]:
    if not PROGRESS_PATH.exists():
        return set()
    try:
        with PROGRESS_PATH.open(encoding="utf-8") as f:
            return set(json.load(f).get("done", []))
    except (json.JSONDecodeError, OSError):
        return set()


def _save_progress(done: set[int]) -> None:
    with PROGRESS_PATH.open("w", encoding="utf-8") as f:
        json.dump({"done": sorted(done), "count": len(done)}, f)


# =============================================================================
# 主流程
# =============================================================================
def _parse_args(argv: list[str]) -> dict[str, Any]:
    opts = {"limit": None, "input": INPUT_PATH, "output": OUTPUT_PATH, "smiles": None}
    i = 0
    while i < len(argv):
        if argv[i] == "--input" and i + 1 < len(argv):
            opts["input"] = ROOT / argv[i + 1]
            i += 2
        elif argv[i] == "--output" and i + 1 < len(argv):
            opts["output"] = ROOT / argv[i + 1]
            i += 2
        elif argv[i] == "--limit" and i + 1 < len(argv):
            opts["limit"] = int(argv[i + 1])
            i += 2
        elif argv[i] == "--smiles" and i + 1 < len(argv):
            opts["smiles"] = argv[i + 1]
            i += 2
        else:
            print(f"[WARN] ignoring unknown arg: {argv[i]}", file=sys.stderr)
            i += 1
    return opts


def _run_single(smiles: str) -> int:
    print(f"SMILES: {smiles}", file=sys.stderr)
    result = _call_llm(f"SMILES: {smiles}")
    if result is None:
        print("ERROR: LLM returned no parseable result", file=sys.stderr)
        return 1
    entry = _sanitize(result, smiles)
    for key in ("id", "error"):
        entry.pop(key, None)
    print(json.dumps(entry, ensure_ascii=False, indent=2))
    return 0


def main() -> None:
    opts = _parse_args(sys.argv[1:])
    if opts["smiles"]:
        sys.exit(_run_single(opts["smiles"]))
    input_path, output_path = opts["input"], opts["output"]
    limit = opts["limit"]

    with input_path.open(encoding="utf-8") as f:
        rows = json.load(f)
    total = len(rows)
    if limit is not None:
        total = min(total, limit)
        rows = rows[:total]
    print(f"Loaded {total} entries from {input_path.name}")

    payload = _load_results()
    done = _load_done()
    print(f"Resuming: {len(done)} already done, {total - len(done)} remaining")

    cfg = _load_llm_config()
    concurrency = max(1, int(cfg.get("concurrency", 1)))
    n_ok, n_err = 0, 0
    batch = 0
    pool = ThreadPoolExecutor(max_workers=concurrency)
    try:
        fut_to_i: dict = {}
        for i, row in enumerate(rows):
            if i in done:
                continue
            fut_to_i[pool.submit(_call_llm, _build_user_message(row), cfg)] = i
        for fut in as_completed(fut_to_i):
            i = fut_to_i[fut]
            row = rows[i]
            smiles = row.get("smiles", "")
            print(f"[{i + 1}/{total}] {smiles[:48]!r} ...", end=" ", flush=True)

            result = fut.result()
            entry = _sanitize(result, smiles)
            entry["id"] = row.get("id")
            payload["results"][str(i)] = entry
            done.add(i)

            if entry["error"]:
                n_err += 1
                print("ERR")
            else:
                n_ok += 1
                print("OK")

            batch += 1
            if batch % SAVE_EVERY == 0:
                payload["meta"] = {"total": total, "processed": len(done)}
                _save_results(payload)
                _save_progress(done)
                print(f"  [saved {len(done)}/{total}]")
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        print("\nInterrupted. Saving progress...")
    finally:
        pool.shutdown(wait=True)

    payload["meta"] = {"total": total, "processed": len(done)}
    _save_results(payload)
    _save_progress(done)
    print(f"\nDone: OK={n_ok} ERR={n_err}, saved {len(done)}/{total} -> {output_path.name}")

    if len(done) >= total:
        PROGRESS_PATH.unlink(missing_ok=True)
        print("All entries processed — progress file removed.")


if __name__ == "__main__":
    main()
