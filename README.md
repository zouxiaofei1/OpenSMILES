# namepredict

Rule-based SMILES → bilingual (English / Chinese) IUPAC namer, with an agent loop for iterative rule refinement and a local ChemAgent Console.

## 安装 / 环境

```bash
python -m venv .venv
# Windows (Git Bash)
source .venv/Scripts/activate
# Linux / macOS
# source .venv/bin/activate

pip install -e ".[dev]"
```

Windows notes (Git Bash): always `source .venv/Scripts/activate` before `python` / `pytest`. Paths in this repo use forward slashes and work under Git Bash.

## 快速开始

```bash
pip install -e ".[dev]"
python tools/merge_datasets.py --tiers smiles_tiers.json --chebi chebi20_test_1k.json --out data/merged_benchmark.json
pytest -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 100
uvicorn server.app:app --host 127.0.0.1 --port 8765
# 循环（无 UI）
python -m agent_loop.loop --max-iters 1 --mock-pi
```

Console URL: [http://127.0.0.1:8765/](http://127.0.0.1:8765/)

### Agent loop CLI

```bash
python -m agent_loop.loop --help
python -m agent_loop.loop --max-iters 1 --mock-pi
```

Optional flags (aligned with `LoopConfig`):

| Flag | Meaning |
|------|---------|
| `--max-iters N` | Max cycles (default 100) |
| `--mock-pi` | Use `MockPiRunner` (no real `pi` CLI) |
| `--bench-limit N` | Limit benchmark rows per cycle |
| `--data PATH` | Benchmark JSON (default `data/merged_benchmark.json`) |
| `--k N` | Stop after N consecutive no-improve cycles |
| `--target-dual F` | Stop when dual accuracy ≥ F |

Without `--mock-pi`, the loop uses the real `PiRunner` (requires a working local `pi`).

**Important:** even with `--mock-pi`, each cycle still runs **real** benchmark + git gate on the current working tree (lint, pytest, commit of allowlisted paths). Revert is **allowlist-only**: it restores `src/namepredict/`, `tests/unit/`, and `skills/chem-tdd-skill/` to the cycle base SHA and cleans untracked files under those prefixes — it does **not** `git reset --hard` the whole repo, so dirty tracked files outside the allowlist survive. Prefer a clean allowlist tree; untracked `data/` / docs are usually fine. Console **Start** shows the same warning before launching.

## Tests

```bash
pytest -q
```

## Benchmark

```bash
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 100
```

If `data/merged_benchmark.json` is missing, generate it with `tools/merge_datasets.py` as in 快速开始.

## Console / API

```bash
uvicorn server.app:app --host 127.0.0.1 --port 8765
```

Open [http://127.0.0.1:8765/](http://127.0.0.1:8765/) for the ChemAgent Console (sessions, namer probe, loop control, SSE logs). API binds to localhost by default.

Namer tab supports drawing structures via self-hosted Ketcher (`web/vendor/ketcher/`, version in `VERSION`). Rebuild vendor:

```bash
bash tools/build_ketcher_vendor.sh
```

Requires Node 18+. Text SMILES naming works even if vendor is missing.

## Package layout

```
src/namepredict/   # installable namer (Layer0–5 + cache)
agent_loop/        # self-improve loop: config, git gate, pi runner, state
server/            # FastAPI app + routes (name, loop, sessions, events)
web/               # ChemAgent Console static UI (dark OLED)
tools/             # merge_datasets, structure_lint, fail_cluster, …
benchmarks/        # bilingual field-wise scorer CLI
skills/            # chem-tdd-skill and related agent skills
tests/             # pytest suite
data/              # merged_benchmark.json (generated) and dataset copies
```

## License / status

Research / internal tooling for SMILES→IUPAC rule development. See `docs/` for architecture and API notes.
