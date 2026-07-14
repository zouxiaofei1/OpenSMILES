"""Scan benchmark fails for fused/heteroarene parent keywords."""
from __future__ import annotations

import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "src"))

from benchmarks.benchmark import _is_fail, _load_rows  # noqa: E402
from benchmarks.benchmark_parallel import (  # noqa: E402
    _chunksize,
    _default_workers,
    _init_worker,
    _score_row,
)


def main() -> None:
    rows = _load_rows(Path("data/merged_benchmark.json"), None)
    workers = _default_workers()
    cs = _chunksize(len(rows), workers)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        results = list(ex.map(_score_row, rows, chunksize=cs))

    fails = []
    for score, pred_en, pred_zh, row in results:
        if _is_fail(score):
            fails.append(
                {
                    "smiles": row.get("smiles"),
                    "en": row.get("english_name") or row.get("en") or "",
                    "zh": row.get("chinese_name") or row.get("zh") or "",
                    "pred_en": pred_en,
                    "pred_zh": pred_zh,
                }
            )
    print("fails", len(fails))
    keys = [
        "phenanthrene",
        "phthalazine",
        "cinnoline",
        "quinazolin",
        "quinoxalin",
        "anthracene",
        "naphthalene",
        "acridine",
        "purine",
        "pteridine",
        "phenazine",
        "phenothiazine",
        "isoquinoline",
        "quinoline",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "coumarin",
        "chromene",
        "chroman",
        "xanthene",
        "fluorene",
        "biphenyl",
        "stilbene",
        "phenanthridine",
        "benzo[",
        "acridon",
        "carbazole",
        "dibenzofuran",
        "dibenzothiophene",
    ]
    buckets: dict[str, list] = defaultdict(list)
    other = []
    for f in fails:
        en = (f["en"] or "").lower()
        hit = next((k for k in keys if k in en), None)
        if hit:
            buckets[hit].append(f)
        elif any(
            x in en
            for x in (
                "benzo",
                "naphth",
                "pheno",
                "quin",
                "indol",
                "azol",
                "azine",
                "furan",
                "thio",
                "imidazo",
                "pyrazolo",
            )
        ):
            other.append(f)

    print("--- by keyword ---")
    for k, v in sorted(buckets.items(), key=lambda kv: -len(kv[1])):
        print(f"{k}: {len(v)}")
        for s in v[:4]:
            print("  gold:", s["en"][:100])
            print("  pred:", (s["pred_en"] or "")[:100])
            print("  smi :", (s["smiles"] or "")[:70])
    print("other", len(other))
    for s in other[:20]:
        print("  gold:", s["en"][:100])
        print("  pred:", (s["pred_en"] or "")[:100], "|", (s["smiles"] or "")[:50])


if __name__ == "__main__":
    main()
