"""Send each benchmark entry to OpenAI LLM for substituent/ring classification.

Usage:
    python scripts/llm_classify_benchmark.py

Output:
    tools/benchmark_llm_classified.jsonl   – one JSON object per line (append-safe)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any
import os
# =============================================================================
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          CONFIGURATION                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
# =============================================================================

# ── OpenAI API settings ──
OPENAI_API_KEY  = "sk-c8af56a13775460ab46c1618c8e9f91c"       # <-- 改成你的 key
OPENAI_BASE_URL = "https://api.deepseek.com"  # <-- 如需代理/中转改这里
OPENAI_MODEL    = "deepseek-v4-pro"                     # <-- 模型名称

# ── Runtime settings ──
LLM_DELAY       = 0.3   # 请求间隔（秒），避免触发速率限制
SAVE_EVERY      = 2    # 每 N 条保存一次进度

# ── Paths ──
ROOT        = Path(__file__).resolve().parents[1]
INPUT_PATH  = ROOT / "scripts" / "benchmark_pred_preview_data.json"
OUTPUT_PATH = ROOT / "scripts" / "benchmark_llm_classified.jsonl"
PROGRESS_PATH = ROOT / "scripts" / "benchmark_llm_progress.json"

# =============================================================================
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          END OF CONFIGURATION                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
# =============================================================================


# ---------------------------------------------------------------------------
# OpenAI client (lazy init)
# ---------------------------------------------------------------------------
_client = None


def _get_client():
    global _client
    if _client is None:
        from openai import OpenAI
        _client = OpenAI(api_key=OPENAI_API_KEY, base_url=OPENAI_BASE_URL)
    return _client


# ---------------------------------------------------------------------------
# Prompt
# ---------------------------------------------------------------------------
# Allowed substituent names — must match exactly one of these keys.
# These are the only substituent types recognised by the downstream naming pipeline.
_ALLOWED_SUBSTITUENTS: set[str] = {
    # ── Simple alkyl (systematic) ──
    "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl", "heptyl", "octyl", "nonyl", "decyl",
    # ── Branched alkyl (retained) ──
    "isopropyl", "isobutyl", "sec-butyl", "tert-butyl", "neopentyl", "isopentyl",
    # ── Alkenyl / alkynyl ──
    "vinyl", "allyl", "isopropenyl", "propargyl", "crotyl", "cinnamyl",
    "ethenyl", "ethynyl", "prop-1-enyl", "prop-2-ynyl",
    # ── Aryl / aralkyl ──
    "phenyl", "benzyl", "tolyl", "phenethyl", "benzhydryl", "trityl",
    # ── Heteroaryl ──
    "furyl", "thienyl", "pyridyl", "quinolyl", "isoquinolyl", "anthryl", "phenanthryl",
    "furfuryl", "thenyl",
    # ── Cycloalkyl ──
    "cyclopropyl", "cyclobutyl", "cyclopentyl", "cyclohexyl", "adamantyl",
    # ── Oxygen-containing ──
    "hydroxy", "methoxy", "ethoxy", "propoxy", "butoxy", "phenoxy",
    "hydroperoxy", "oxo", "formyloxy", "acetoxy",
    # ── Sulfur-containing ──
    "sulfanyl", "methylsulfanyl", "ethylsulfanyl",
    "methylsulfinyl", "ethylsulfinyl", "methylsulfonyl", "ethylsulfonyl",
    "sulfo", "mesyl", "tosyl", "triflyl",
    "chlorosulfonyl", "isothiocyanato",
    # ── Nitrogen-containing ──
    "amino", "methylamino", "dimethylamino", "ethylamino", "diethylamino",
    "nitroso", "nitro", "azido", "hydrazinyl", "diazenyl", "diazo",
    "anilino", "cyano", "isocyano",
    "carbamoyl", "sulfamoyl", "amidino", "guanidino", "ureido",
    "hydrazono", "hydroxyimino", "imino",
    # ── Halogen / haloalkyl ──
    "fluoro", "chloro", "bromo", "iodo",
    "trifluoromethyl", "trichloromethyl", "tribromomethyl",
    "difluoromethyl", "pentafluoroethyl",
    # ── P / B / Se ──
    "phosphanyl", "phosphono", "phosphonooxy",
    "boranyl", "borono", "selanyl",
    # ── Acyl ──
    "formyl", "acetyl", "benzoyl", "propionyl", "butyryl",
    "isobutyryl", "valeryl", "oxamoyl",
    "methoxycarbonyl", "ethoxycarbonyl",
    "chlorocarbonyl",
    # ── Other heteroatom ──
    "silyl", "trimethylsilyl", "germyl", "stannyl",
}

_ALLOWED_NAMES_STR = ", ".join(sorted(_ALLOWED_SUBSTITUENTS))

SYSTEM_PROMPT = f"""\
You are a cheminformatics assistant. Given a molecule's SMILES and its IUPAC names \
(English and Chinese), identify:

1. substituents: a JSON object mapping each substituent type name to its count.
   **CRITICAL**: every substituent name MUST be chosen from this exact list:
   {_ALLOWED_NAMES_STR}

   If a substituent you observe is not in this list, find the closest standard name \
   from the list (e.g. "isopentyl" not "3-methylbutyl"; "methylsulfanyl" not "methylthio"; \
   "oxo" not "keto"; "methoxycarbonyl" not "carbomethoxy").
   If there are no substituents, return an empty object {{}}.

2. ring_count: total number of rings in the molecule (0 for acyclic).

Return ONLY a valid JSON object, no extra text or markdown fences.
Example: {{"substituents": {{"methyl": 2, "hydroxy": 1}}, "ring_count": 1}}\
"""


def _build_user_message(row: dict[str, Any]) -> str:
    parts: list[str] = []
    if row.get("s"):
        parts.append(f"SMILES: {row['s']}")
    if row.get("ge"):
        parts.append(f"English IUPAC name (gold): {row['ge']}")
    if row.get("gz"):
        parts.append(f"Chinese IUPAC name (gold): {row['gz']}")
    if row.get("en"):
        parts.append(f"English predicted name: {row['en']}")
    if row.get("zh"):
        parts.append(f"Chinese predicted name: {row['zh']}")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# LLM call
# ---------------------------------------------------------------------------
def _normalize_substituents(raw: dict[str, Any]) -> dict[str, int]:
    """Filter and normalise substituent dict: keep only allowed names, coerce counts to int."""
    out: dict[str, int] = {}
    for k, v in raw.items():
        k = k.strip().lower()
        if k in _ALLOWED_SUBSTITUENTS:
            try:
                out[k] = int(v)
            except (ValueError, TypeError):
                out[k] = 1
        else:
            # Try fuzzy: case-insensitive match
            match = next((a for a in _ALLOWED_SUBSTITUENTS if a.lower() == k), None)
            if match:
                try:
                    out[match] = int(v)
                except (ValueError, TypeError):
                    out[match] = 1
            else:
                print(f"  [WARN] Unknown substituent '{k}' discarded", file=sys.stderr)
    return out


def _call_llm(user_message: str, model: str) -> dict[str, Any] | None:
    client = _get_client()
    try:
        resp = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": user_message},
            ],
            temperature=0.0,
            max_tokens=512,
            response_format={"type": "json_object"},
        )
        raw = resp.choices[0].message.content.strip()
        return json.loads(raw)
    except json.JSONDecodeError:
        # Try to extract JSON from markdown fences
        import re
        m = re.search(r"\{[^{}]*\"substituents\"[^{}]*\}", raw, re.DOTALL)
        if m:
            try:
                return json.loads(m.group(0))
            except json.JSONDecodeError:
                pass
        print(f"  [WARN] Failed to parse JSON response: {raw[:200]}", file=sys.stderr)
        return None
    except Exception as exc:
        print(f"  [ERR] API call failed: {exc}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Progress tracking (append-safe resume)
# ---------------------------------------------------------------------------
def _load_progress() -> set[int]:
    """Return set of already-processed indices."""
    if not PROGRESS_PATH.exists():
        return set()
    try:
        with PROGRESS_PATH.open(encoding="utf-8") as f:
            data = json.load(f)
        return set(data.get("done_indices", []))
    except (json.JSONDecodeError, OSError):
        return set()


def _save_progress(done_indices: set[int]) -> None:
    with PROGRESS_PATH.open("w", encoding="utf-8") as f:
        json.dump({"done_indices": sorted(done_indices), "count": len(done_indices)}, f)


def _append_result(entry: dict[str, Any]) -> None:
    with OUTPUT_PATH.open("a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    model = os.environ.get("OPENAI_MODEL", "deepseek-v4-pro")
    delay = float(os.environ.get("LLM_DELAY", "0.3"))  # seconds between requests

    # Load input
    print(f"Loading: {INPUT_PATH}")
    with INPUT_PATH.open(encoding="utf-8") as f:
        rows = json.load(f)
    total = len(rows)
    print(f"Total entries: {total}")

    # Resume
    done = _load_progress()
    if done:
        print(f"Resuming: {len(done)} already done, {total - len(done)} remaining")

    # Process
    n_success, n_fail = 0, 0
    try:
        for i, row in enumerate(rows):
            if i in done:
                continue

            smiles = row.get("s", "")[:60]
            print(f"[{i + 1}/{total}] {smiles} ...", end=" ", flush=True)

            user_msg = _build_user_message(row)
            result = _call_llm(user_msg, model)

            if result is None:
                n_fail += 1
                print("FAIL")
                # Still record a failure entry so the index aligns
                _append_result({
                    "idx": i,
                    "smiles": row.get("s"),
                    "error": True,
                    "substituents": {},
                    "ring_count": None,
                })
            else:
                n_success += 1
                entry = {
                    "idx": i,
                    "smiles": row.get("s"),
                    "substituents": _normalize_substituents(result.get("substituents", {})),
                    "ring_count": result.get("ring_count"),
                }
                _append_result(entry)
                print("OK")

            done.add(i)

            # Save progress periodically
            if (i + 1) % SAVE_EVERY == 0:
                _save_progress(done)

            if LLM_DELAY > 0:
                time.sleep(LLM_DELAY)

    except KeyboardInterrupt:
        print("\nInterrupted. Saving progress...")
        _save_progress(done)
        print(f"Progress saved. {len(done)}/{total} done. Resume by re-running the script.")
        sys.exit(0)
    finally:
        _save_progress(done)

    # Summary
    print(f"\nDone! Success: {n_success}, Failed: {n_fail}")
    print(f"Output: {OUTPUT_PATH}")

    # Clean up progress file on full completion
    if len(done) >= total:
        PROGRESS_PATH.unlink(missing_ok=True)
        print("All entries processed — progress file removed.")


if __name__ == "__main__":
    main()
