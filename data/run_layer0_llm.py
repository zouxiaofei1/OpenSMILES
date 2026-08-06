"""LLM-driven Layer0 分析器：对 merged_benchmark.json 每个测试例调用 OpenAI-compatible API，
让 LLM 产出 layer0 层的结果（SMILES 解析 / 无机-有机分类 / 盐解离元数据）。

行为：
  - 每轮处理一个测试例，串行直到全部完成；
  - 每处理 SAVE_EVERY 条（默认 5）把完整结果写回 layer0_output.json；
  - 进度记录在独立的状态文件，中断（Ctrl-C）后重跑即从断点继续；
  - 全部完成后保留结果文件、删除进度文件。

用法（在项目根目录）:
    ./.venv/Scripts/python.exe data/run_layer0_llm.py
    ./.venv/Scripts/python.exe data/run_layer0_llm.py --limit 50     # 只跑前 50 条（测试）
    ./.venv/Scripts/python.exe data/run_layer0_llm.py --input data/merged_benchmark_shuffled.json
    ./.venv/Scripts/python.exe data/run_layer0_llm.py --smiles 'CC(=O)O[Na]'   # 单物质：只打印结果，不写文件

环境变量（覆盖默认配置）:
    OPENAI_API_KEY   默认 sk-c8af56a13775460ab46c1618c8e9f91c
    OPENAI_BASE_URL  默认 https://api.deepseek.com
    OPENAI_MODEL     默认 deepseek-v4-pro
    LLM_DELAY        默认 0.3（秒，两次请求间隔，避免限流）
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
SAVE_EVERY      = 5          # 每处理多少条写盘一次
TIMEOUT         = 60.0       # 单次请求超时（秒）
MAX_TOKENS      = 4000       # deepseek-v4-pro 为 reasoning 模型，token 同时覆盖思考+回答，需给足
CALL_RETRIES    = 3          # 空响应/解析失败时的重试次数

ROOT        = Path(__file__).resolve().parents[1]
SETTINGS_PATH = ROOT / "server" / "settings.json"
INPUT_PATH  = ROOT / "data" / "merged_benchmark.json"
OUTPUT_PATH = ROOT / "data" / "layer0_output.json"
PROGRESS_PATH = ROOT / "data" / "layer0_progress.json"


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

# =============================================================================
# LLM 客户端（httpx 直连 /chat/completions，无需安装 openai 包）
# =============================================================================
def _chat_completions(messages: list[dict], *, model: str, api_key: str, base_url: str) -> str:
    """POST /chat/completions，返回 assistant 的文本内容。非 2xx 抛异常。"""
    payload = {
        "model": model,
        "messages": messages,
        "temperature": 0.0,
        "max_tokens": MAX_TOKENS,
    }
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    url = base_url.rstrip("/") + "/chat/completions"
    # 多数兼容端点支持 response_format；不支持时报错后去掉重试一次
    for with_format in (True, False):
        p = dict(payload)
        if with_format:
            p["response_format"] = {"type": "json_object"}
        with httpx.Client(timeout=TIMEOUT) as client:
            resp = client.post(url, headers=headers, json=p)
            if resp.status_code == 200:
                data = resp.json()
                return data["choices"][0]["message"]["content"]
            if not with_format or "response_format" not in resp.text.lower():
                resp.raise_for_status()
    raise RuntimeError("unreachable")

# =============================================================================
# PROMPT
# =============================================================================
SYSTEM_PROMPT = """\
You are the Layer-0 preprocessor of a SMILES→IUPAC naming pipeline. For the given SMILES \
and (optional) gold IUPAC names, perform exactly three Layer-0 duties:

1. parseable — can the SMILES be parsed into a valid molecule? (true/false)
2. classification — is it organic, inorganic, or ambiguous?
   - inorganic: no carbon skeleton / common inorganic species (H2O, CO2, salts of \
     inorganic acids, metal oxides/halides, noble gases, etc.)
3. salt — if it is a salt or an acid/amine addduct, record it:
   - alkali metal salt (Li/Na/K cation): {"type": "alkali_metal", "metal": "sodium", \
     "metal_zh": "钠", "n_metal": 1}
   - hydrochloride (HCl adduct or [Cl-] counterion): {"type": "hydrochloride", \
     "acid": "hydrochloride", "acid_zh": "盐酸盐"}
   - any other salt/adduct (ammonium, Ca/Mg/Fe, sulfate/nitrate/phosphate, hydrate, \
     co-crystal...): {"type": "other", "description": "short description"}
   - if not a salt: null

Also output:
- organic_smiles: the organic fragment's SMILES after stripping salt counterions \
  (keep the original string mostly intact; null if purely inorganic or unparseable).
- n_fragments: number of independent molecular units (dot-separated; salt with \
  counterion counts as >=2).
- notes: one short sentence justifying the classification.

Respond ONLY with a single JSON object, no markdown fences, no extra text. Example:
{"parseable": true, "classification": "organic", "salt": {"type": "alkali_metal", \
"metal": "sodium", "metal_zh": "钠", "n_metal": 1}, "organic_smiles": "CC(=O)[O-]", \
"n_fragments": 2, "notes": "Sodium acetate: Na+ counterion detached."}\
"""


def _build_user_message(row: dict[str, Any]) -> str:
    parts = [f"SMILES: {row.get('smiles', '')}"]
    if row.get("english_name"):
        parts.append(f"English IUPAC name (gold): {row['english_name']}")
    if row.get("chinese_name"):
        parts.append(f"Chinese IUPAC name (gold): {row['chinese_name']}")
    return "\n".join(parts)


def _call_llm(user_message: str, cfg: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """调用 LLM 并解析 JSON；空响应/解析失败重试 CALL_RETRIES 次。

    cfg 为 None 时每次调用现读 server/settings.json（保证设置页改动即时生效）。
    reasoning 模型的思考过程会占用 max_tokens，思考一长 content 就可能被截断
    为空或残缺，所以不能一次失败就放弃。
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


# =============================================================================
# 字段规整：保证输出结构完整、类型正确
# =============================================================================
def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return default


def _sanitize(raw: dict[str, Any] | None, smiles: str) -> dict[str, Any]:
    raw = raw or {}
    classification = str(raw.get("classification", "ambiguous")).strip().lower()
    if classification not in ("organic", "inorganic", "ambiguous"):
        classification = "ambiguous"
    salt = raw.get("salt")
    if isinstance(salt, dict):
        salt_type = str(salt.get("type", "other"))
        if salt_type not in ("alkali_metal", "hydrochloride", "other"):
            salt_type = "other"
        salt = {"type": salt_type, **salt}
    return {
        "id": None,
        "smiles": smiles,
        "parseable": _as_bool(raw.get("parseable"), True),
        "classification": classification,
        "salt": salt,
        "organic_smiles": raw.get("organic_smiles") if raw.get("organic_smiles") else None,
        "n_fragments": int(raw.get("n_fragments", 1) or 1),
        "notes": str(raw.get("notes", "") or ""),
        "error": not isinstance(raw, dict),
    }


# =============================================================================
# 状态与结果存取
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
# MAIN
# =============================================================================
def _parse_args(argv: list[str]) -> dict[str, Any]:
    opts = {"limit": None, "input": INPUT_PATH, "output": OUTPUT_PATH, "progress": PROGRESS_PATH, "smiles": None}
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
    """单物质模式：调 LLM 做 layer0 分析，只把结果打印到控制台，不写文件。"""
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
    input_path, output_path, progress_path = opts["input"], opts["output"], opts["progress"]
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
        progress_path.unlink(missing_ok=True)
        print("All entries processed — progress file removed.")


if __name__ == "__main__":
    main()
