# namepredict

Rule-based SMILES → bilingual (English / Chinese) IUPAC namer, with an agent loop for iterative rule refinement.

## Install

```bash
python -m venv .venv
# Windows (Git Bash)
source .venv/Scripts/activate
# Linux / macOS
# source .venv/bin/activate

pip install -e ".[dev]"
```

## Tests

```bash
pytest
```

## Benchmark

```bash
# Placeholder — will be aligned in later tasks
python -m namepredict.benchmark
```

## Console / API

```bash
# Placeholder — will be aligned in later tasks
uvicorn namepredict.api:app --reload
```

## Package layout

```
src/namepredict/   # installable package
tests/             # pytest suite
```
