# namepredict Agent 总提示（守则 + 工作流）

> **读者**：实现 Agent（pi / 人工会话）、编排器作者、本机 Console 操作者。  
> **仓库根**：`E:\dev\chem`（Windows；路径用 `pathlib`，勿写死 Unix-only）。  
> **目标**：基于 **规则**（非 ML）把 SMILES 命名为 **英文 + 中文 IUPAC**；用 **可门禁的自循环** 迭代规则，直至 dual 达标或停机。

每轮动态上下文（iter / dual / 失败簇 / bench 命令）由编排器注入 `agent_loop/prompts/cycle.md`。  
**本文件是静态总则**；与 `cycle.md`、`skills/chem-tdd-skill/` 冲突时，以 **更严者** 为准。

---

## 1. 使命与成功标准

| 项 | 说明 |
|----|------|
| 产品 | `namepredict`：`SMILESNNamer().name(smiles) → NameResult(en, zh, success, …)` |
| 方法 | **仅规则 + 受限 cache**；禁止用深度学习或外部 LLM **生成化学名** |
| 数据 | 主集 `data/merged_benchmark.json`（tiers + chebi 合并） |
| 计分 | **分字段**：有 gold 才计；英文轻度归一化；中文 strip 后严格相等；**dual** = 应测语言全对 |
| 循环目标 | dual 可配置（默认 **≥ 99%**），或用户中断 / 连续 K=5 无提升 / max-iters |
| 提升定义 | dual↑，**或** dual 相同且 fails↓ |
| 回退 | 门禁失败或未提升 → 编排器 **仅恢复 allowlist 路径**（不做整库 `reset --hard`） |

**非目标（Agent 勿做）**：补全 chebi 中文金标、改计分公式、公网多用户鉴权、10ms 性能专项（除非人类单独立项）。

---

## 2. 系统架构（你在哪一层）

```
ChemAgent Console (web/)  127.0.0.1:8765
        │ HTTP + SSE
        ▼
FastAPI (server/)  loop / sessions / name / events / settings / bench
        │
        ▼
agent_loop (Python 编排器)
  状态机：bench → cluster → 渲染 cycle 提示 → pi → lint → pytest → bench → git 门禁 → session
        │ 每轮 prompt + cwd
        ▼
pi CLI（真循环；需 ANTHROPIC_API_KEY）
        │ 只改 allowlist 内代码
        ▼
src/namepredict/  Layer0–5 + 入口 + cache
```

| 组件 | 路径 | 职责 |
|------|------|------|
| 命名管道 | `src/namepredict/` | 生产命名逻辑 |
| 计分 | `benchmarks/benchmark.py` | **只读金标**；分字段计分 |
| 合并数据 | `tools/merge_datasets.py` | 生成 merged 主集 |
| 结构 lint | `tools/structure_lint.py` | 行数、函数体、层 import、禁 SMILES 特判、cache≤100 |
| 失败簇 | `tools/fail_cluster.py` | 按 feature 聚类，建议 Layer |
| 编排 | `agent_loop/` | 循环、git 门禁、pi、状态、secrets |
| TDD Skill | `skills/chem-tdd-skill/` | 每轮强制 TDD |
| 控制台 | `server/` + `web/` | 启停、Namer、Bench、API Key |

**实现 Agent 只负责「改命名规则 + 单测」**；commit/revert、启停、全量 bench 调度由 **编排器 / 人类** 负责。

---

## 3. 命名管道（Layer0–5）

入口：`SMILESNNamer.name` → cache → preprocess → analyze → select_parent → extract_substituents → number → assemble。

| Layer | 目录 | 职责 |
|-------|------|------|
| L0 | `layer0/` | 预处理：解析 SMILES、sanitize、盐/单元拆分等 |
| L1 | `layer1/` | 结构分析：环、链、官能团等 |
| L2 | `layer2/` | 母体选择（IUPAC 优先级 / 保留名等） |
| L3 | `layer3/` | 取代基抽取与类型 |
| L4 | `layer4/` | 编号与位次 |
| L5 | `layer5/` | 中英文名组装 |
| cache | `cache/` | **仅**俗名/例外，**≤ 100** 条 |

- **禁止**在一个 layer 包内混入另一 layer 的职责文件。  
- **禁止**层向上错误依赖（由 `structure_lint` 检查）。  
- 规则依据：`docs/iupac/`、`docs/cleaned/`、`docs/iupac/cn/`（只读参考）。

---

## 4. 硬性守则（违反 → 本轮作废 / 门禁失败）

### 4.1 方法与数据

1. **禁止 ML/LLM 命名**：不得用模型直接输出 IUPAC 名；LLM 只允许改 **代码/测试**（由 pi 会话完成）。  
2. **cache ≤ 100**：俗名+例外；禁止把 benchmark 金标整表搬进 cache。  
3. **禁止 SMILES 特判**：`if smiles == "..."`、完整 SMILES→名称大哈希表等单分子捷径。  
4. **禁止改金标与计分**：  
   - 勿改：`data/*`、`chebi20_test_1k.json`、`smiles_tiers.json`、`data/merged_benchmark.json`  
   - 勿改：`benchmarks/benchmark.py` 的计分逻辑  
5. **禁止**为过门禁而删测试、弱化断言、放宽归一化、改 lint 放行名单（除非人类明确要求）。

### 4.2 代码形态

6. **可写路径（生产）**：仅  
   - `src/namepredict/`（入口、`layer0`–`layer5`、受限 `cache`）  
   - `tests/unit/`  
   - 编排器要求时的 `agent_loop/memory/` 小结（一般由编排器写 session）  
7. **只读**：`benchmarks/`、`tools/`（可 **运行** 不可改逻辑）、`docs/` 规则原文、`server/`、`web/`、`agent_loop` 编排核心（除非人类任务明确改编排）。  
8. **单文件 ≤ 500 行**；**单函数/方法体 ≤ 10 行**（可拆私有辅助函数）。  
9. **Windows**：`pathlib`；命令在 Git Bash 下；`source .venv/Scripts/activate`。

### 4.3 粒度 S3（N=3）

10. **默认**每轮 **1 个主规则**。  
11. **上限**：主规则 1 + 同簇同层附属 ≤ **N−1**，合计规则点 **≤ 3**。  
12. 禁止「顺手修一片」；其余失败留给后续轮。  
13. 连续同分且同主簇 ≥ 2 轮：优先 **L1/L2/L4** 结构修复，**禁止**只在 L5 打补丁。

### 4.4 Git 与安全

14. **Agent 不要自己 `git commit` / `git reset`**；由 `GitGate` 决定 commit 或 allowlist-only 恢复。  
15. **不要**把 API Key 写入代码、测试、session、progress、日志。  
16. Key 仅存在本机 `agent_loop/memory/secrets.json`（gitignore）；Console 配置后注入 pi 子进程 `ANTHROPIC_API_KEY`。

---

## 5. 整个工作流

### 5.1 一图总览

```
[人类/Console]
  配 API Key → Start 循环（或 CLI）
       │
       ▼
┌──────────────── agent_loop 一轮 ────────────────┐
│ 1. load STATE.json                              │
│ 2. baseline bench（可 limit；子进程冷导入）        │
│ 3. cluster_failures → 主簇摘要                   │
│ 4. 渲染 cycle.md → .cycle_current.md            │
│ 5. pi -p  ← 注入本文件精神 + skill + 簇上下文    │
│ 6. structure_lint(src/namepredict)              │
│ 7. pytest tests/unit                            │
│ 8. bench again                                  │
│ 9. gate: 提升? → commit allowlist : revert 路径  │
│10. 写 session JSON + progress.md                │
│11. 停机? STOP / no_improve≥K / dual≥目标 / max  │
└─────────────────────────────────────────────────┘
```

### 5.2 编排器状态机（`agent_loop/loop.py`）

1. **load state** — `agent_loop/memory/STATE.json`（iter, best_dual, last_dual, no_improve, …）  
2. **baseline bench** — 子进程 `python -m benchmarks.benchmark --json`（避免旧 bytecode）  
3. **cluster** — `tools.fail_cluster.cluster_failures`  
4. **render prompt** — `agent_loop/prompts/cycle.md` 变量：  
   `{{iter}} {{dual}} {{cluster}} {{skill_path}} {{constraints}} {{bench_cmd}}`  
5. **pi run** — `PiRunner`：`pi -p`，prompt 走 stdin；env 可含 `ANTHROPIC_API_KEY`  
6. **structure_lint**  
7. **pytest**  
8. **bench**  
9. **gate**  
   - lint/pytest 失败 → revert  
   - 未提升 → revert + `no_improve++`  
   - 提升且有 allowlist 脏文件 → commit  
   - 提升但无脏文件 → `noop`（仍计 no_improve）  
   - **last_dual（能力 dual）**：revert/noop 时保持 dual_before，避免假达标停机  
10. **session** — `agent_loop/memory/sessions/cycle-XXXX.json`  
11. **stop** — `STOP` 文件 / `no_improve ≥ K` / `last_dual ≥ target` / `iter ≥ max_iters`

### 5.3 实现 Agent 单轮（你在 pi 里要做的）

**严格顺序（TDD）** — 详见 `skills/chem-tdd-skill/SKILL.md`：

1. 读失败簇、建议 Layer、`docs/cleaned` / `docs/iupac`、可选 `progress.md`  
2. 抽 **3–10 正例 + ≥1 近邻负例**  
3. **先写** `tests/unit/` 测试；文件头：  
   `# IUPAC: P-…`  
   `# Layer: Lx`  
4. **pytest 必须先红**（未实现时失败；禁止假绿）  
5. **再改** Layer/入口/cache 使测试绿  
6. pytest 绿 → `python tools/structure_lint.py --root src/namepredict` → 跑编排器给出的 `{{bench_cmd}}`  
7. 写 **回合小结**（主规则、文件、红绿证据、lint/bench 摘要、未覆盖项）  
8. **不要 git commit**

### 5.4 双语断言约定

```python
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

r = SMILESNNamer().name(smiles)
assert r.success
assert normalize_en(r.en) == normalize_en(en_gold)
# 仅当有中文金标时：
assert normalize_zh(r.zh) == normalize_zh(zh_gold)
```

- 英文：strip / 小写 / 空白合并（以 `normalize_en` 为准）  
- 中文：strip 后严格相等  
- 无中文金标：跳过中文断言，勿瞎编 zh

### 5.5 本机 Console / CLI（操作者）

```bash
source .venv/Scripts/activate
pip install -e ".[dev]"

# 数据（若缺 merged）
python tools/merge_datasets.py \
  --tiers smiles_tiers.json --chebi chebi20_test_1k.json \
  --out data/merged_benchmark.json

pytest -q
python tools/structure_lint.py --root src/namepredict
python -m benchmarks.benchmark --data data/merged_benchmark.json --limit 50

# 控制台
uvicorn server.app:app --host 127.0.0.1 --port 8765
# 浏览器 http://127.0.0.1:8765/
# 顶栏设置 Anthropic API Key → 才能 Start（真 pi）
# Bench 页可手动全量测试（阶段日志 SSE）

# 无 UI 循环
python -m agent_loop.loop --max-iters 1 --mock-pi --bench-limit 30
# 真 pi：不要 --mock-pi，且本机/Console 已配置 Key
```

| 操作 | 说明 |
|------|------|
| API Key | 存 `agent_loop/memory/secrets.json`；未配置 → Start 400 / 弹窗 |
| Start | 真 pi + 真实 lint/pytest/bench/git 门禁；确认后运行 |
| 手动 Bench | 全量 `merged_benchmark.json`；与 loop 可并行，占 CPU |
| Namer | 本地规则，**不需要** Key |

---

## 6. 常用命令速查

| 目的 | 命令 |
|------|------|
| 单测 | `pytest -q` / `pytest tests/unit -q` |
| 结构 | `python tools/structure_lint.py --root src/namepredict` |
| 计分 | `python -m benchmarks.benchmark --data data/merged_benchmark.json [--limit N] [--json]` |
| 循环 | `python -m agent_loop.loop --help` |
| 服务 | `uvicorn server.app:app --host 127.0.0.1 --port 8765` |

---

## 7. 记忆与产物（只读/可写边界）

| 路径 | 谁写 | 内容 |
|------|------|------|
| `agent_loop/memory/STATE.json` | 编排器 | 迭代状态 |
| `agent_loop/memory/progress.md` | 编排器 / 约定格式日志 | 轮次摘要 |
| `agent_loop/memory/sessions/` | 编排器 | 每轮 JSON |
| `agent_loop/memory/STOP` | Console/人类 | 停机 |
| `agent_loop/memory/secrets.json` | Console/Settings | API Key（**gitignore**） |
| `agent_loop/prompts/.cycle_current.md` | 编排器 | 本轮完整提示 |

Agent **可读** progress / session / 簇摘要；**勿**把密钥写进上述文件。

---

## 8. 回合小结模板（Agent 输出）

1. **主规则**（IUPAC ID）与 **Layer**  
2. **改动文件**（相对仓库根）  
3. **测试**：正例/负例数量；文件名  
4. **红→绿证据**（实现前失败、实现后通过）  
5. **structure_lint / pytest / benchmark** 摘要  
6. **未覆盖失败**（留给后续）  
7. **风险 / follow-up**

---

## 9. 相关文件索引

| 文档 | 用途 |
|------|------|
| `prompt.md`（本文件） | 总则：守则 + 全工作流 |
| `agent_loop/prompts/cycle.md` | 每轮动态注入模板 |
| `skills/chem-tdd-skill/SKILL.md` | 强制 TDD 细节 |
| `skills/chem-tdd-skill/checklist.md` | 检查清单 |
| `docs/superpowers/specs/2026-07-12-smiles-iupac-agent-loop-design.md` | 系统设计 |
| `docs/superpowers/plans/2026-07-12-smiles-iupac-agent-loop.md` | 脚手架实现计划（已完成） |
| `docs/superpowers/specs/2026-07-12-console-apikey-bench-design.md` | Console Key + 全量 Bench |
| `README.md` | 人类快速开始 |

---

## 10. 一句话纪律

**先测后码、先红后绿；规则进 Layer，不进金标；一小步 S3，门禁说了算；密钥不落库，commit 不越权。**
