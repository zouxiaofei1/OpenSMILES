# SMILES→IUPAC Agent Loop + ChemAgent Console Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建成基于规则的 SMILES→中英文 IUPAC 命名骨架、合并 benchmark 门禁、TDD skill、pi 自循环编排器，以及 Dark OLED 的 ChemAgent Console 控制台。

**Architecture:** Python 包 `namepredict`（Layer0–5 流水线）产出 `{en,zh}`；`benchmarks` + `tools` 负责计分与结构 lint；`agent_loop` 以失败簇驱动循环（TDD→改代码→lint→bench→git commit/reset）；FastAPI `server` + `web` 提供本机控制面（session 侧栏、启停、Namer、SSE 日志）。设计规格见 `docs/superpowers/specs/2026-07-12-smiles-iupac-agent-loop-design.md`。

**Tech Stack:** Python 3.11+、RDKit、pytest、FastAPI、uvicorn、SSE；前端静态 HTML + Tailwind CDN + Lucide；可选本机 `pi` CLI；git 门禁。

## Global Constraints

- 禁止 ML 命名；俗名/例外缓存 **≤100** 条
- 生产代码仅 `src/namepredict/` 入口 + `layer0`–`layer5` + 受限 `cache`
- 单文件 **≤500** 行；单函数/方法 **≤10** 行
- 禁止 `if smiles == "..."` 单题特判；禁止改 benchmark 金标与计分逻辑
- 英文比较：轻度归一化；中文：strip 后严格相等；分字段计分
- 粒度 S3：默认 1 主规则，同簇同层最多 N=3
- 分数下降 **>0.5%** → `git reset --hard`；提升 = dual↑ 或 同分 fails↓
- 停机：用户中断 / 连续 K=5 无提升 / dual≥99% / max-iters
- Console 默认绑定 **127.0.0.1**
- UI：Dark OLED，无 emoji 图标，对比度 ≥4.5:1，focus 可见
- 所有用户可见文案与提交说明可用中文；代码标识符英文
- 工作目录：`E:\dev\chem`（Windows；路径用 pathlib）

---

## File Structure (lock-in)

```
E:\dev\chem\
├── pyproject.toml
├── requirements.txt
├── README.md
├── data\
│   ├── smiles_tiers.json          # copy or link from root originals
│   ├── chebi20_test_1k.json
│   └── merged_benchmark.json      # generated
├── src\namepredict\
│   ├── __init__.py
│   ├── namer.py
│   ├── constants.py
│   ├── types.py
│   ├── layer0\preprocessor.py
│   ├── layer1\analyzer.py
│   ├── layer2\parent_selector.py
│   ├── layer3\substituent_extractor.py
│   ├── layer4\numbering.py
│   ├── layer5\assembler.py
│   └── cache\common_names.py
├── tests\
│   ├── conftest.py
│   └── unit\
│       ├── test_namer_smoke.py
│       ├── test_normalize.py
│       └── test_cache_limit.py
├── benchmarks\
│   ├── __init__.py
│   └── benchmark.py
├── tools\
│   ├── merge_datasets.py
│   ├── structure_lint.py
│   └── fail_cluster.py
├── agent_loop\
│   ├── __init__.py
│   ├── config.py
│   ├── state.py
│   ├── events.py
│   ├── git_gate.py
│   ├── pi_runner.py
│   ├── loop.py
│   ├── prompts\cycle.md
│   └── memory\  (progress.md, STATE.json, sessions\)
├── skills\chem-tdd-skill\
│   ├── SKILL.md
│   ├── checklist.md
│   └── templates\test_rule.py.tmpl
├── server\
│   ├── __init__.py
│   ├── app.py
│   ├── deps.py
│   ├── routes_loop.py
│   ├── routes_sessions.py
│   ├── routes_name.py
│   └── routes_events.py
└── web\
    ├── index.html
    ├── css\app.css
    └── js\app.js
```

---

### Task 1: 仓库初始化与 Python 包装

**Files:**
- Create: `pyproject.toml`
- Create: `requirements.txt`
- Create: `README.md`
- Create: `src/namepredict/__init__.py`
- Create: `.gitignore`

**Interfaces:**
- Produces: installable package name `namepredict` (src layout); pytest discovers `tests/`

- [ ] **Step 1: 初始化 git（若尚未）并写 `.gitignore`**

```gitignore
__pycache__/
*.py[cod]
.venv/
venv/
.env
.pytest_cache/
.mypy_cache/
*.egg-info/
dist/
build/
data/merged_benchmark.json
agent_loop/memory/sessions/
agent_loop/prompts/.cycle_current.md
agent_loop/memory/STATE.json
*.log
.DS_Store
Thumbs.db
pi/node_modules/
```

Run:
```bash
cd /e/dev/chem
git status || git init
```

- [ ] **Step 2: 写 `pyproject.toml` 与 `requirements.txt`**

`pyproject.toml`:
```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "namepredict"
version = "0.1.0"
description = "Rule-based SMILES to bilingual IUPAC namer"
readme = "README.md"
requires-python = ">=3.11"
dependencies = [
  "rdkit>=2023.9.1",
  "fastapi>=0.110.0",
  "uvicorn[standard]>=0.27.0",
  "pydantic>=2.0",
  "httpx>=0.27.0",
]

[project.optional-dependencies]
dev = ["pytest>=8.0", "pytest-cov>=4.0"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
```

`requirements.txt`:
```
rdkit>=2023.9.1
fastapi>=0.110.0
uvicorn[standard]>=0.27.0
pydantic>=2.0
httpx>=0.27.0
pytest>=8.0
```

- [ ] **Step 3: 最小包入口与 README**

`src/namepredict/__init__.py`:
```python
"""Rule-based SMILES → bilingual IUPAC naming."""

__version__ = "0.1.0"
```

`README.md`: 简述目标、安装 `pip install -e ".[dev]"`、跑 `pytest`、跑 benchmark、起 Console 的命令（可先写占位命令，后续任务对齐）。

- [ ] **Step 4: 创建 venv 并安装**

```bash
cd /e/dev/chem
python -m venv .venv
source .venv/Scripts/activate  # Git Bash on Windows
pip install -e ".[dev]"
python -c "import namepredict; print(namepredict.__version__)"
```
Expected: `0.1.0`

- [ ] **Step 5: Commit**

```bash
git add .gitignore pyproject.toml requirements.txt README.md src/namepredict/__init__.py
git commit -m "chore: init namepredict package and gitignore"
```

---

### Task 2: 核心类型、归一化与缓存上限

**Files:**
- Create: `src/namepredict/types.py`
- Create: `src/namepredict/constants.py`
- Create: `src/namepredict/cache/common_names.py`
- Create: `src/namepredict/cache/__init__.py`
- Create: `tests/unit/test_normalize.py`
- Create: `tests/unit/test_cache_limit.py`
- Create: `tests/conftest.py`

**Interfaces:**
- Produces:
  - `normalize_en(name: str) -> str`
  - `normalize_zh(name: str) -> str`
  - `NameResult` dataclass: `en, zh, success, source, time_ms, meta`
  - `CommonNameCache` with `get(smiles) -> NameResult | None`, `max_entries=100`

- [ ] **Step 1: 写失败测试 `tests/unit/test_normalize.py`**

```python
from namepredict.constants import normalize_en, normalize_zh

def test_normalize_en_lower_and_space():
    assert normalize_en("  Ethanol  ") == "ethanol"
    assert normalize_en("propan-2-one") == "propan-2-one"
    assert normalize_en("A  B") == "a b"

def test_normalize_zh_strip_only():
    assert normalize_zh("  乙醇  ") == "乙醇"
    assert normalize_zh("乙醇") == "乙醇"
```

- [ ] **Step 2: 跑测确认失败**

```bash
pytest tests/unit/test_normalize.py -v
```
Expected: FAIL import error

- [ ] **Step 3: 实现 `constants.py` 与 `types.py`**

`src/namepredict/types.py`:
```python
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

@dataclass
class NameResult:
    en: str
    zh: str
    success: bool
    source: str = "iupac"
    time_ms: float = 0.0
    meta: dict[str, Any] = field(default_factory=dict)
```

`src/namepredict/constants.py`:
```python
from __future__ import annotations
import re

_WS = re.compile(r"\s+")

def normalize_en(name: str) -> str:
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",").replace(" ,", ",")
    return s

def normalize_zh(name: str) -> str:
    return (name or "").strip()
```

- [ ] **Step 4: 缓存测试与实现**

`tests/unit/test_cache_limit.py`:
```python
import pytest
from namepredict.cache.common_names import CommonNameCache
from namepredict.types import NameResult

def test_cache_rejects_over_100():
    c = CommonNameCache(max_entries=100)
    for i in range(100):
        c.put(f"C{i}", NameResult(en=f"n{i}", zh=f"名{i}", success=True, source="cache"))
    with pytest.raises(ValueError):
        c.put("overflow", NameResult(en="x", zh="x", success=True, source="cache"))
```

`src/namepredict/cache/common_names.py`:
```python
from __future__ import annotations
from namepredict.types import NameResult

class CommonNameCache:
    def __init__(self, max_entries: int = 100) -> None:
        self.max_entries = max_entries
        self._data: dict[str, NameResult] = {}

    def get(self, smiles: str) -> NameResult | None:
        return self._data.get(smiles)

    def put(self, smiles: str, result: NameResult) -> None:
        if smiles in self._data:
            self._data[smiles] = result
            return
        if len(self._data) >= self.max_entries:
            raise ValueError(f"cache exceeds max_entries={self.max_entries}")
        self._data[smiles] = result

    def __len__(self) -> int:
        return len(self._data)
```

- [ ] **Step 5: 跑通测试并提交**

```bash
pytest tests/unit/test_normalize.py tests/unit/test_cache_limit.py -v
git add src/namepredict tests
git commit -m "feat: add normalize helpers, NameResult, cache limit"
```

---

### Task 3: Layer0–5 骨架 + Namer 管道（最小规则）

**Files:**
- Create: each `layerN/*.py` and `__init__.py`
- Create: `src/namepredict/namer.py`
- Create: `tests/unit/test_namer_smoke.py`

**Interfaces:**
- Consumes: `NameResult`, `CommonNameCache`, RDKit
- Produces: `SMILESNNamer.name(smiles: str) -> NameResult`
- Layer functions (pure-ish, each public function ≤10 lines; split helpers as needed):
  - `layer0.preprocess(smiles) -> Mol | None`
  - `layer1.analyze(mol) -> dict`
  - `layer2.select_parent(info) -> dict`
  - `layer3.extract_substituents(info, parent) -> list`
  - `layer4.number(parent, substituents) -> dict`
  - `layer5.assemble(...) -> NameResult`

**最小可命名集合（规则，非特判 SMILES 字符串表）:**
- 直链烷烃 C1–C10 母体名表（按碳数）
- 一元醇：最长碳链 + 羟基后缀 → ethanol/乙醇 等（至少覆盖 `C`, `CC`, `CCO`）

- [ ] **Step 1: 失败测试**

```python
# tests/unit/test_namer_smoke.py
from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh

def test_methane():
    r = SMILESNNamer().name("C")
    assert r.success
    assert normalize_en(r.en) == "methane"
    assert normalize_zh(r.zh) == "甲烷"

def test_ethanol():
    r = SMILESNNamer().name("CCO")
    assert r.success
    assert normalize_en(r.en) == "ethanol"
    assert normalize_zh(r.zh) == "乙醇"
```

- [ ] **Step 2: 跑测确认失败**

```bash
pytest tests/unit/test_namer_smoke.py -v
```

- [ ] **Step 3: 实现各层（保持函数短小）**

实现要点（写入代码时拆成 ≤10 行函数）：

`layer0/preprocessor.py`: `Chem.MolFromSmiles` → `AddHs` 可选 → 失败返回 None  

`layer1/analyzer.py`: 原子数、键、是否含羟基（氧连 H 与 C）、碳链长度启发式  

`layer2/parent_selector.py`: 有羟基则母体为连羟基的碳链；否则最长碳链  

`layer3/substituent_extractor.py`: 骨架阶段返回空列表  

`layer4/numbering.py`: 一元醇羟基位尽量为 1（乙醇可省略位次，按 PIN 习惯输出 ethanol）  

`layer5/assembler.py`:  
- 烷烃英文 `methane..decane` / 中文 `甲..癸 + 烷`  
- 醇：`stem + "ol"` / `干 + "醇"`（ethanol/乙醇 特例词干表按碳数，不是 SMILES 特判）

`namer.py`:
```python
def name(self, smiles: str) -> NameResult:
    # 1 cache 2 preprocess 3 analyze 4 parent 5 subst 6 number 7 assemble
    # 计时 time_ms；失败 success=False en/zh 空或 unknown
```

- [ ] **Step 4: 跑测通过**

```bash
pytest tests/unit/test_namer_smoke.py -v
```
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/namepredict tests/unit/test_namer_smoke.py
git commit -m "feat: layer0-5 pipeline with alkane/alcohol smoke names"
```

---

### Task 4: 合并数据集

**Files:**
- Create: `tools/merge_datasets.py`
- Create: `data/` copies or path refs to root JSON
- Generate: `data/merged_benchmark.json`
- Create: `tests/unit/test_merge_datasets.py`（可用临时小 fixture，不依赖全量）

**Interfaces:**
- Produces CLI: `python tools/merge_datasets.py --tiers PATH --chebi PATH --out PATH`
- Record fields: `id, smiles, english_name, chinese_name, source, tier, features, eval_en, eval_zh`

- [ ] **Step 1: 写 merge 单测（小 fixture）**

```python
# tests/unit/test_merge_datasets.py
from pathlib import Path
import json
from tools.merge_datasets import merge_records

def test_merge_prefers_zh_and_sets_flags(tmp_path: Path):
    tiers = [{"id": 1, "smiles": "CCO", "english_name": "ethanol", "chinese_name": "乙醇", "tier": 1, "features": ["alcohol"]}]
    chebi = [{"id": 9, "smiles": "CCO", "english_name": "ethanol", "chinese_name": "", "tier": 1, "features": []}]
    out = merge_records(tiers, chebi)
    assert len(out) == 1
    assert out[0]["eval_zh"] is True
    assert out[0]["chinese_name"] == "乙醇"
    assert out[0]["source"] == "smiles_tiers"
```

注意：若 `tools` 非包，在 `merge_datasets.py` 提供可 import 的 `merge_records`，测试里把项目根加入 path（conftest）或改成 `from importlib`）。

- [ ] **Step 2: 实现 merge**

逻辑：
1. 读 JSON 列表  
2. RDKit canonical SMILES 作键；失败则用原 smiles  
3. 冲突：优先 `chinese_name` 非空；再优先 `source==smiles_tiers`  
4. `eval_en = bool(english_name.strip())`；`eval_zh = bool(chinese_name.strip())`  
5. id：`tiers-{id}` / `chebi-{id}`

- [ ] **Step 3: 跑全量合并**

```bash
mkdir -p data
cp smiles_tiers.json data/ 2>/dev/null || true
cp chebi20_test_1k.json data/ 2>/dev/null || true
python tools/merge_datasets.py \
  --tiers /e/dev/chem/smiles_tiers.json \
  --chebi /e/dev/chem/chebi20_test_1k.json \
  --out /e/dev/chem/data/merged_benchmark.json
python -c "import json;d=json.load(open('data/merged_benchmark.json',encoding='utf-8'));print(len(d), sum(x['eval_zh'] for x in d), sum(x['eval_en'] for x in d))"
```
Expected: 条数 ≤ 810+3297；`eval_zh` 约 770 左右；`eval_en` 接近全量

- [ ] **Step 4: Commit 脚本与（可选）不提交巨型 merged 若 gitignore——保留生成命令**

```bash
git add tools/merge_datasets.py tests/unit/test_merge_datasets.py
git commit -m "feat: merge tiers+chebi into benchmark dataset"
```

---

### Task 5: Benchmark 计分脚本

**Files:**
- Create: `benchmarks/__init__.py`
- Create: `benchmarks/benchmark.py`
- Create: `tests/unit/test_benchmark_score.py`

**Interfaces:**
- Produces:
  - `score_record(pred_en, pred_zh, row) -> dict`
  - `run_benchmark(namer, data_path, limit=None) -> Report`
  - CLI: `python -m benchmarks.benchmark --data data/merged_benchmark.json [--limit 50]`
- Report keys: `acc_en, acc_zh, acc_dual, n_en, n_zh, n_dual, fails, by_source, by_tier`

- [ ] **Step 1: 计分单测**

```python
from benchmarks.benchmark import score_record

def test_dual_only_en_when_no_zh():
    row = {"english_name": "ethanol", "chinese_name": "", "eval_en": True, "eval_zh": False}
    r = score_record("ethanol", "", row)
    assert r["en_ok"] is True
    assert r["zh_ok"] is None
    assert r["dual_ok"] is True

def test_dual_requires_both_when_flags():
    row = {"english_name": "ethanol", "chinese_name": "乙醇", "eval_en": True, "eval_zh": True}
    r = score_record("ethanol", "酒精", row)
    assert r["en_ok"] is True
    assert r["zh_ok"] is False
    assert r["dual_ok"] is False
```

- [ ] **Step 2: 实现 `benchmarks/benchmark.py`**

使用 `normalize_en` / `normalize_zh`；调用 `SMILESNNamer().name`；打印百分比：

```
en=..% (a/b) zh=..% (c/d) dual=..% (e/f) fails=N
```

禁止修改输入 JSON。

- [ ] **Step 3: 跑小样本**

```bash
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 20
```
Expected: 无崩溃；甲烷/乙醇若在前部可能命中

- [ ] **Step 4: Commit**

```bash
git add benchmarks tests/unit/test_benchmark_score.py
git commit -m "feat: bilingual field-wise benchmark scorer"
```

---

### Task 6: structure_lint

**Files:**
- Create: `tools/structure_lint.py`
- Create: `tests/unit/test_structure_lint.py`

**Interfaces:**
- CLI exit 0=pass, 1=fail
- Checks under `src/namepredict/`:
  1. 文件行数 ≤500  
  2. 函数体行数 ≤10（用 ast 统计 FunctionDef/AsyncFunctionDef 体非空行）  
  3. 缓存 dict 字面量或 `CommonNameCache` 使用处不强制；另检 `common_names` 模块内预置条目 ≤100  
  4. 禁止匹配 `if\s+smiles\s*==`  
  5. layer 目录不得 import 更高层（简单：layer0 不得 import layer1–5 等）

- [ ] **Step 1: 单测用临时目录构造违规文件**

- [ ] **Step 2: 实现 lint 并在当前骨架上跑通**

```bash
python tools/structure_lint.py --root src/namepredict
```
Expected: exit 0（若函数超长则拆分后重跑）

- [ ] **Step 3: Commit**

```bash
git add tools/structure_lint.py tests/unit/test_structure_lint.py
git commit -m "feat: structure lint for layers, line limits, no smiles special-case"
```

---

### Task 7: 失败簇工具

**Files:**
- Create: `tools/fail_cluster.py`
- Create: `tests/unit/test_fail_cluster.py`

**Interfaces:**
- `cluster_failures(fails: list[dict], top_k=5) -> list[Cluster]`
- Cluster: `key, size, features, sample_smiles, suggest_layers, error_langs`

- [ ] **Step 1: 测试聚类按 feature 聚合**

- [ ] **Step 2: 实现：优先 features 元组；空 features 用 `source`+`error_langs`**

- [ ] **Step 3: Commit**

```bash
git add tools/fail_cluster.py tests/unit/test_fail_cluster.py
git commit -m "feat: failure clustering for agent loop"
```

---

### Task 8: chem-tdd-skill

**Files:**
- Create: `skills/chem-tdd-skill/SKILL.md`
- Create: `skills/chem-tdd-skill/checklist.md`
- Create: `skills/chem-tdd-skill/templates/test_rule.py.tmpl`

**Interfaces:**
- Produces: Agent 可读 skill；编排器每轮注入路径

- [ ] **Step 1: 写 `SKILL.md`（完整条款，中文）**

必须包含：
1. 先测后码；先红后绿  
2. 从失败簇抽 3–10 SMILES + ≥1 负例  
3. 双语断言规则  
4. 文件头 `IUPAC:` / `Layer:`  
5. S3 边界 N=3  
6. 禁止单 SMILES 特判与改金标  
7. 命令：`pytest tests/unit -q`、`python tools/structure_lint.py`、`python -m benchmarks.benchmark ...`  
8. 回合小结格式

- [ ] **Step 2: checklist + 测试模板**

模板示例：
```python
# IUPAC: P-XX.X
# Layer: Lx
import pytest
from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh

CASES = [
    # ("smiles", "en", "zh_or_None"),
]

@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_rule(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
```

- [ ] **Step 3: Commit**

```bash
git add skills/chem-tdd-skill
git commit -m "docs: add chem-tdd-skill for agent cycles"
```

---

### Task 9: agent_loop 配置、状态、事件、git 门禁

**Files:**
- Create: `agent_loop/config.py`
- Create: `agent_loop/state.py`
- Create: `agent_loop/events.py`
- Create: `agent_loop/git_gate.py`
- Create: `tests/unit/test_git_gate.py`
- Create: `tests/unit/test_state.py`

**Interfaces:**
- `LoopConfig`: `k=5, target_dual=0.99, max_iters=100, drop_tol=0.005, n_rules=3, data_path, mock_pi, cwd`
- `StateStore`: load/save `STATE.json`；`progress.md` append
- `EventBus`: `subscribe(fn)`, `publish(type, payload)`
- `GitGate`: `snapshot() -> sha`, `commit(msg, paths)`, `revert(sha)`, `allowed_paths`

- [ ] **Step 1: state/event 单测（tmp_path）**

- [ ] **Step 2: git_gate 单测** — 在 tmp git 仓 touch 文件 commit/revert

实现 `git_gate.py` 用 subprocess `git`；Windows 兼容。

- [ ] **Step 3: Commit**

```bash
git add agent_loop tests/unit/test_git_gate.py tests/unit/test_state.py
git commit -m "feat: loop config, state store, event bus, git gate"
```

---

### Task 10: pi_runner + cycle prompt 模板

**Files:**
- Create: `agent_loop/pi_runner.py`
- Create: `agent_loop/prompts/cycle.md`
- Create: `tests/unit/test_pi_runner.py`

**Interfaces:**
- `PiRunner.run(prompt_path: Path, cwd: Path, timeout_s: int) -> PiResult(success, log, exit_code)`
- `MockPiRunner`：不调 pi，写一个 no-op 或按环境变量失败/成功
- Prompt 模板变量：`{{iter}} {{dual}} {{cluster}} {{skill_path}} {{constraints}} {{bench_cmd}}`

- [ ] **Step 1: Mock runner 测试**

- [ ] **Step 2: 实现真实 runner**

尝试：
```bash
pi -p --print < prompt_file
```
若 CLI 不同，读 `pi/docs/usage.md` 校正；失败时 `PiResult(success=False)`。

- [ ] **Step 3: 写 `cycle.md` 完整提示词（中文）**

嵌入：硬约束、S3、TDD 流程、只改允许路径、完成后自列小结。

- [ ] **Step 4: Commit**

```bash
git add agent_loop/pi_runner.py agent_loop/prompts/cycle.md tests/unit/test_pi_runner.py
git commit -m "feat: pi runner and cycle prompt template"
```

---

### Task 11: 主循环 `loop.py`

**Files:**
- Create: `agent_loop/loop.py`
- Create: `tests/unit/test_loop_gate.py`
- Create: `agent_loop/memory/progress.md`（初始空结构）

**Interfaces:**
- `AgentLoop(config, namer_factory, pi_runner, bus).run_once() -> CycleResult`
- `AgentLoop.run() -> None` 直到停机
- `CycleResult`: iter, dual_before, dual_after, decision, session_path

**状态机（必须实现）:**
1. load state  
2. baseline bench（可 limit 配置加速测试）  
3. cluster failures  
4. render prompt  
5. pi run  
6. structure_lint  
7. pytest  
8. bench  
9. gate → commit or reset  
10. write session json + progress  
11. stop checks：STOP 文件、`state.no_improve >= K`、`dual >= target`、`iter >= max`

- [ ] **Step 1: 门禁单测（mock pi + mock bench 分数）**

```python
def test_gate_reverts_on_dual_drop(tmp_path):
    # before dual=0.5 after=0.4 (>0.5% drop) => decision revert
    ...

def test_gate_commits_on_dual_up(tmp_path):
    ...
```

- [ ] **Step 2: 实现 `loop.py`（拆短函数）**

- [ ] **Step 3: dry-run 一轮**

```bash
python -c "from agent_loop.loop import AgentLoop; from agent_loop.config import LoopConfig; from agent_loop.pi_runner import MockPiRunner; AgentLoop(LoopConfig(max_iters=1, mock_pi=True, bench_limit=30), pi_runner=MockPiRunner()).run()"
```
Expected: 写出 session；不崩溃；可能 revert（mock 未改代码）

- [ ] **Step 4: Commit**

```bash
git add agent_loop tests/unit/test_loop_gate.py
git commit -m "feat: agent loop state machine with bench git gate"
```

---

### Task 12: FastAPI Server

**Files:**
- Create: `server/app.py`
- Create: `server/deps.py`
- Create: `server/routes_loop.py`
- Create: `server/routes_sessions.py`
- Create: `server/routes_name.py`
- Create: `server/routes_events.py`
- Create: `tests/unit/test_api_name.py`

**Interfaces:**
- `GET /api/v1/loop/state`
- `POST /api/v1/loop/start|pause|stop`
- `GET /api/v1/sessions` / `GET /api/v1/sessions/{iter}`
- `POST /api/v1/name` body `{"smiles":"CCO"}`
- `GET /api/v1/events` SSE
- Host: `127.0.0.1:8765`

- [ ] **Step 1: `test_api_name` 用 httpx ASGITransport**

```python
from fastapi.testclient import TestClient
from server.app import app

def test_name_ethanol():
    c = TestClient(app)
    r = c.post("/api/v1/name", json={"smiles": "CCO"})
    assert r.status_code == 200
    body = r.json()
    assert body["success"] is True
    assert "ethanol" in body["en"].lower()
```

- [ ] **Step 2: 实现路由**

`loop` 控制：模块级 `LoopController`（线程跑 `AgentLoop.run`）；`pause` 设标志；`stop` 写 `agent_loop/memory/STOP` 并可选 `force`。

SSE：从 `EventBus` 队列读事件，`text/event-stream`。

- [ ] **Step 3: 手动起服务**

```bash
uvicorn server.app:app --host 127.0.0.1 --port 8765
curl -s http://127.0.0.1:8765/api/v1/loop/state
curl -s -X POST http://127.0.0.1:8765/api/v1/name -H "content-type: application/json" -d "{\"smiles\":\"C\"}"
```

- [ ] **Step 4: Commit**

```bash
git add server tests/unit/test_api_name.py
git commit -m "feat: FastAPI control plane for loop, sessions, name, SSE"
```

---

### Task 13: ChemAgent Console 前端

**Files:**
- Create: `web/index.html`
- Create: `web/css/app.css`
- Create: `web/js/app.js`

**Interfaces:**
- Consumes: Task 12 APIs
- UI: 按 `design-system/chemagent-console/MASTER.md` + `pages/console.md`

- [ ] **Step 1: 搭建 Dark OLED 壳**

- CSS 变量：
  - `--color-background: #0F172A`
  - `--color-foreground: #F8FAFC`
  - `--color-accent: #22C55E`
  - `--color-destructive: #EF4444`
  - `--color-muted: #272F42`
  - `--color-border: #475569`
- 字体 Inter（Google fonts + `font-display: swap`）
- 布局：顶栏 + 左 280px session + 主 tabs
- Lucide CDN 或 inline SVG（Play / Pause / Square / Flask）

- [ ] **Step 2: 实现 JS**

- `fetchState()` 轮询 2s 或 SSE 优先  
- Session 列表渲染；点击加载 `/sessions/{iter}` 到 Cycle/Prompt  
- Start/Pause/Stop 按钮；Stop 用 `confirm`  
- Namer 表单 POST `/name`  
- Logs：SSE `log.line` append；`prefers-reduced-motion` 时关闭脉冲动画  

- [ ] **Step 3: 静态挂载**

在 `server/app.py`：
```python
app.mount("/", StaticFiles(directory="web", html=True), name="web")
```
注意：API 路由注册在 mount 之前。

- [ ] **Step 4: 浏览器手测清单**

- 打开 `http://127.0.0.1:8765/`  
- Namer 输入 `CCO` 得到 ethanol/乙醇  
- 侧栏空状态文案可见  
- 键盘 Tab 有焦点环  
- 顶栏 dual/状态显示  

- [ ] **Step 5: Commit**

```bash
git add web server/app.py
git commit -m "feat: ChemAgent Console dark OLED UI with namer and sessions"
```

---

### Task 14: 端到端接线与文档

**Files:**
- Modify: `README.md`
- Modify: `agent_loop/config.py` defaults if needed
- Create: `docs/superpowers/plans/2026-07-12-smiles-iupac-agent-loop.md`（本文件已存在则更新验收节）

- [ ] **Step 1: README 命令一节**

```markdown
## 快速开始
pip install -e ".[dev]"
python tools/merge_datasets.py --tiers smiles_tiers.json --chebi chebi20_test_1k.json --out data/merged_benchmark.json
pytest -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 100
uvicorn server.app:app --host 127.0.0.1 --port 8765
# 循环（无 UI）
python -m agent_loop.loop --max-iters 1 --mock-pi
```

- [ ] **Step 2: 全量回归**

```bash
pytest -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 50
```

- [ ] **Step 3: 最终 commit**

```bash
git add README.md
git commit -m "docs: wire up quickstart for namer, loop, and console"
```

---

## Spec Coverage Checklist

| Spec 项 | Task |
|---------|------|
| 合并数据集 + 分字段计分 | 4, 5 |
| 英归一化 / 中严格 | 2, 5 |
| Layer0–5 + 入口 | 3 |
| 硬约束 lint | 6 |
| 失败簇驱动 | 7, 11 |
| TDD skill S3 | 8, 11 |
| git commit/revert 门禁 | 9, 11 |
| pi runner + mock | 10, 11 |
| 停机条件 | 11 |
| session 快照 | 11, 12 |
| Console UI + API | 12, 13 |
| 设计系统 dark OLED | 13 |
| 缓存 ≤100 | 2, 6 |

## Placeholder Scan

无 TBD/TODO 步骤；CLI 若与本机 pi 不符，在 Task 10 按 `pi/docs/usage.md` 校正具体参数（仍须实现 Mock 路径）。

## Type Consistency

- `NameResult` 字段贯穿 namer / cache / API  
- `acc_dual` 为 0–1 float；UI 显示百分比  
- Session 文件字段与 spec §11 一致  

---

**Plan complete.** 保存于 `docs/superpowers/plans/2026-07-12-smiles-iupac-agent-loop.md`。
