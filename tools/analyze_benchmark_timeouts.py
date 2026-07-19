from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "src"))


def _run_one(smiles: str, conn) -> None:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    from namepredict.namer import SMILESNNamer

    mol = Chem.MolFromSmiles(smiles)
    features = {
        "heavy_atoms": mol.GetNumHeavyAtoms() if mol else None,
        "bonds": mol.GetNumBonds() if mol else None,
        "rings": mol.GetRingInfo().NumRings() if mol else None,
        "fragments": len(Chem.GetMolFrags(mol)) if mol else None,
    }
    namer = SMILESNNamer()
    conn.send(("ready", features))
    started = time.perf_counter()
    try:
        result = namer.name(smiles)
        conn.send(("result", {
            "elapsed_s": time.perf_counter() - started,
            "success": result.success,
            "en": result.en,
            "zh": result.zh,
            "meta": result.meta,
        }))
    except Exception as exc:
        conn.send(("error", {
            "elapsed_s": time.perf_counter() - started,
            "error": f"{type(exc).__name__}: {exc}",
        }))
    finally:
        conn.close()


def _measure(index: int, row: dict, timeout_s: float) -> dict:
    smiles = str(row.get("smiles") or "")
    parent_conn, child_conn = mp.Pipe(duplex=False)
    process = mp.Process(target=_run_one, args=(smiles, child_conn))
    process.start()
    child_conn.close()
    features = {}
    try:
        if not parent_conn.poll(30):
            process.terminate()
            process.join()
            return {"index": index, "smiles": smiles, "status": "startup_timeout"}
        kind, features = parent_conn.recv()
        if kind != "ready":
            raise RuntimeError(f"unexpected worker message: {kind}")
        if not parent_conn.poll(timeout_s):
            process.terminate()
            process.join()
            return {
                "index": index,
                "smiles": smiles,
                "status": "timeout",
                "timeout_s": timeout_s,
                **features,
            }
        kind, payload = parent_conn.recv()
        process.join(timeout=1)
        return {
            "index": index,
            "smiles": smiles,
            "status": kind,
            **features,
            **payload,
        }
    except (EOFError, OSError, RuntimeError) as exc:
        if process.is_alive():
            process.terminate()
        process.join()
        return {
            "index": index,
            "smiles": smiles,
            "status": "worker_error",
            "error": f"{type(exc).__name__}: {exc}",
            **features,
        }
    finally:
        parent_conn.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--start", type=int, default=2200)
    parser.add_argument("--stop", type=int, default=2500)
    parser.add_argument("--timeout", type=float, default=1.0)
    parser.add_argument("--workers", type=int, default=19)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = json.loads(Path(args.data).read_text(encoding="utf-8"))
    selected = list(enumerate(rows[args.start:args.stop], args.start))
    results = []
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_measure, index, row, args.timeout): index
            for index, row in selected
        }
        for done, future in enumerate(as_completed(futures), 1):
            results.append(future.result())
            if done % 25 == 0 or done == len(selected):
                print(f"progress {done}/{len(selected)}", flush=True)

    results.sort(key=lambda item: item["index"])
    report = {
        "config": {
            "data": args.data,
            "start": args.start,
            "stop": args.stop,
            "timeout_s": args.timeout,
            "workers": args.workers,
        },
        "wall_s": time.perf_counter() - started,
        "results": results,
    }
    Path(args.output).write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    mp.freeze_support()
    main()
