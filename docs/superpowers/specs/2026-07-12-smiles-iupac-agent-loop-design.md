# SMILES → 中英文 IUPAC 规则转换器 + pi Agent 自循环工作流

**日期**: 2026-07-12  
**状态**: 待用户确认  
**目标**: 基于 Python 规则（非 ML）实现 SMILES→IUPAC（英/中）命名；用外层编排器 + pi agent 无人值守迭代，直至达标或用户中断；并提供 **ChemAgent Console** 网页控制台（session 上下文、启停、命名试跑）。

---

## 1. 背景与目标

### 1.1 问题

- 需要**高准确率**（目标双重准确率 ≥ 99%）的 SMILES → 英文/中文 IUPAC 命名。
- 必须**基于规则**实现 IUPAC 蓝皮书逻辑，禁止深度学习；俗名/例外缓存 ≤ 100。
- 规则极多，单次人工/单会话无法完成 → 需要 **可门禁的自循环 Agent 工作流**。

### 1.2 仓库现状（设计时）

| 资源 | 状态 |
|------|------|
| `docs/iupac`, `docs/cleaned`, `docs/iupac/cn` | 规则原文/清洗稿可用 |
| `smiles_tiers.json` (810) | 分层、中英基本齐全 |
| `chebi20_test_1k.json` (3297) | 复杂分子，中文全空 |
| `prompt.txt` | 旧 v4 约束与 Layer 架构可参考 |
| `pi/` | pi coding agent 源码/SDK |
| 转换器源码 / benchmark 脚本 | **缺失，需新建** |

### 1.3 成功标准

1. `name("CCO")` → `en=ethanol`, `zh=乙醇`（及更广规则覆盖）。
2. 合并主集上：**双重准确率 ≥ 99%**（可配置），或用户中断 / 连续 K 轮无提升停机。
3. 自循环：分析失败 → TDD → 改代码 → lint → bench → commit/回退 → 写记忆 → 下一轮。
4. 结构硬约束始终可被自动检查。
5. ChemAgent Console：可查看每轮 session 上下文、启停循环、试跑命名。

---

## 2. 已确认决策

| # | 决策 | 选择 |
|---|------|------|
| 1 | 交付顺序 | 基础设施 → 工作流 → 自动改进 |
| 2 | 数据集 | **合并** tiers + chebi 为新主集 |
| 3 | 计分 | **分字段**：有 gold 才计；双重 = 应测语言全中 |
| 4 | 匹配 | 英文**轻度归一化**；中文**严格相等** |
| 5 | 编排 | **Python 编排器** + 调用 pi CLI/RPC |
| 6 | 循环策略 | **方案 2**：失败簇驱动 + 架构护栏；tier 仅观测 |
| 7 | 代码约束 | **全盘沿用** prompt 硬约束 |
| 8 | 停机 | 用户中断 / 连续 K 轮无提升 / 双重准确率 ≥ 目标 |
| 9 | 粒度 | **S3**：默认 1 主规则；同簇同层最多 N=3 |
| 10 | TDD | benchmark 之外增加 TDD 红→绿门禁 |
| 11 | Skill | **P1**：仓库内 `skills/chem-tdd-skill/`，每轮注入 |
| 12 | 控制台 | **ChemAgent Console** 网页：session 侧栏、启停、Namer 试跑、实时日志 |

**默认参数**: `K=5`，目标双重准确率 `99%`，`max-iters=100`（可配），分数下降 `>0.5%` 回退。  
**提升定义**: 双重准确率上升，**或** 同分但失败条数减少。  
**依赖**: Python 3.11+、RDKit、pytest；包名 `namepredict`。

---

## 3. 架构总览

```
┌──────────────────────────────────────────────────────────────────┐
│  ChemAgent Console (web/)                                        │
│  Session 侧栏 · 启停 · Namer · Bench · Logs · 实时事件             │
└────────────────────────────┬─────────────────────────────────────┘
                             │ HTTP + SSE/WebSocket
                             ▼
                    ┌─────────────────────────────────────┐
                    │  API (FastAPI) + Event Bus          │
                    │  /api/v1/loop|sessions|name|events  │
                    └────────────────┬────────────────────┘
                                     │
                    ┌────────────────▼────────────────────┐
                    │         agent_loop (Python)         │
                    │  状态机 / git 门禁 / 记忆 / 调 pi    │
                    └───────────────┬─────────────────────┘
                                    │ 每轮 prompt + cwd
                                    ▼
                    ┌─────────────────────────────────────┐
                    │              pi agent                 │
                    │  读 skill + 规则 + 失败簇 → 改代码    │
                    └───────────────┬─────────────────────┘
                                    │
          ┌─────────────────────────┼─────────────────────────┐
          ▼                         ▼                         ▼
   structure_lint            pytest (TDD)              benchmark.py
   (层/行数/缓存)             tests/unit                 data/merged
          │                         │                         │
          └─────────────────────────┴─────────────────────────┘
                                    │
                          pass → git commit
                          fail → git reset --hard
```

### 3.1 命名流水线（运行时）

```
SMILES
  → Layer0 Molecule Preprocessor
  → Layer1 Structural Analyzer
  → Layer2 Parent Selector
  → Layer3 Substituent Extractor
  → Layer4 Numbering Engine
  → Layer5 Name Assembler
  → { en, zh, meta }
```

---

## 4. 目录结构（目标）

```
E:\dev\chem\
├── data\
│   ├── smiles_tiers.json          # 只读原始
│   ├── chebi20_test_1k.json       # 只读原始（可从根目录迁移或符号引用）
│   └── merged_benchmark.json      # 合并主集（生成）
├── src\namepredict\
│   ├── __init__.py
│   ├── namer.py                   # 唯一对外入口（可改）
│   ├── constants.py
│   ├── layer0\ ...
│   ├── layer1\ ...
│   ├── layer2\ ...
│   ├── layer3\ ...
│   ├── layer4\ ...
│   ├── layer5\ ...
│   └── cache\                     # 俗名/例外，≤100
├── tests\
│   ├── unit\                      # TDD 单测（按规则/层组织）
│   └── conftest.py
├── benchmarks\
│   └── benchmark.py               # 只读金标；禁止 Agent 改数据
├── tools\
│   ├── merge_datasets.py
│   ├── structure_lint.py
│   └── fail_cluster.py            # 失败簇摘要
├── agent_loop\
│   ├── loop.py                    # 主状态机
│   ├── config.py
│   ├── git_gate.py
│   ├── pi_runner.py               # 调 pi
│   ├── events.py                  # 事件总线 → Console
│   ├── prompts\
│   │   └── cycle.md
│   └── memory\
│       ├── progress.md            # 开发日志 / 完成记录
│       ├── sessions\              # 每轮上下文快照（供侧栏）
│       │   └── cycle-0001.json
│       └── STATE.json             # 轮次、分数、连续无提升计数
├── web\                           # ChemAgent Console 前端
│   ├── index.html
│   ├── assets\
│   └── README.md
├── server\
│   ├── app.py                     # FastAPI 入口
│   ├── routes_loop.py
│   ├── routes_sessions.py
│   ├── routes_name.py
│   └── routes_events.py           # SSE/WebSocket
├── skills\
│   └── chem-tdd-skill\
│       ├── SKILL.md
│       ├── checklist.md
│       └── templates\
│           └── test_rule.py.tmpl
├── design-system\
│   └── chemagent-console\         # ui-ux-pro-max 产出
│       ├── MASTER.md
│       └── pages\console.md
├── docs\                          # 已有 iupac / cleaned / 架构
├── prompt.txt                     # 历史约束参考
├── pyproject.toml / requirements.txt
└── README.md
```

**Agent 可写范围**: `src/namepredict/`（入口 + layer0–5 + cache）、`tests/unit/`、`agent_loop/memory/`（由编排器写为主）。  
**只读**: `data/*` 金标、`benchmarks/benchmark.py`、`docs/iupac/**`、合并脚本生成的主集。

---

## 5. 硬约束（自动 lint）

1. 生产代码仅位于 `src/namepredict/` 的入口与 `layer0`–`layer5`（及受限 `cache`）。
2. **不得**在 layer 文件夹内混入其他 layer 职责；禁止非法跨层 import（允许：高层调用低层公开 API，细节在 lint 白名单中定义）。
3. 单文件 ≤ **500** 行；单函数/方法 ≤ **10** 行。
4. 函数式偏好：输入输出尽量单一；复用已有方法。
5. 禁止「居所取代型」投机捷径；禁止 `if smiles == "..."` 式单题特判。
6. 缓存（俗名+例外）条目 ≤ **100**。
7. 禁止修改 benchmark 金标与 `benchmarks/benchmark.py` 的计分逻辑（工具自身 bugfix 仅人类可做）。
8. 不得使用 ML 模型做命名。

---

## 6. 合并数据集

### 6.1 字段

```json
{
  "id": "tiers-1",
  "smiles": "CCO",
  "english_name": "ethanol",
  "chinese_name": "乙醇",
  "source": "smiles_tiers",
  "tier": 1,
  "features": ["alcohol"],
  "eval_en": true,
  "eval_zh": true
}
```

- 来自 chebi：`id=chebi-{id}`，`eval_zh=false`（中文空），`eval_en` 在有英文名时为 true。
- `source` / `tier` / `features` 保留用于失败聚类与报告，**不**强制路线图。

### 6.2 去重

- 以规范化 SMILES（RDKit canonical）为主键；冲突时优先保留 **含中文 gold** 的条目，其次 `smiles_tiers`。

---

## 7. 计分与归一化

### 7.1 英文轻度归一化（顺序）

1. `strip` + 小写  
2. 合并连续空白为单空格  
3. 统一连字符变体（可选：en-dash → `-`）  
4. 去掉多余空格（如 ` ,` → `,`）  
**不做**：语义等价重排、忽略 locant、同义词映射（同义词仅能通过 ≤100 缓存命中生产侧）。

### 7.2 中文

- `strip` 后**严格相等**（不转繁简、不忽略标点差异）。

### 7.3 指标

- `acc_en` = 正确英文 / `eval_en` 条数  
- `acc_zh` = 正确中文 / `eval_zh` 条数  
- `acc_dual` = 双重命中条数 / 至少测一种语言的条数  
  - 单条双重命中：所有 `eval_*=true` 的语言均正确  

报告同时输出 per-`source`、per-`tier` 分解（观测用）。

---

## 8. 循环状态机

```
INIT
  → LOAD_STATE (progress.md, STATE.json)
  → BENCH_BASELINE
  → while not stop:
        SELECT_CLUSTER      # 失败簇，不按 tier 强制
        BUILD_PROMPT        # 注入 skill + 规则指针 + 约束 + N≤3
        RUN_PI_AGENT        # 期望：TDD 红 → 实现绿 → 自检
        LINT
        PYTEST
        BENCH
        GATE:
          if lint&pytest ok and dual 未降>0.5% and (dual↑ or fails↓):
              COMMIT + LOG + reset no_improve
          else:
              RESET_HARD + LOG fail + no_improve++
        STOP if: user interrupt | no_improve≥K | dual≥target | iter≥max
  → FINAL_REPORT
```

### 8.1 失败簇（SELECT）

输入：错误列表 `(smiles, pred_en, gold_en, pred_zh, gold_zh, features, source)`  
聚类键优先级：`features` 交集 → 错误类型（en/zh/both）→ 简单启发式（含环/含羰基等）  
输出：本轮主簇摘要 + 建议 Layer + 相关 `docs/cleaned` / `docs/iupac` 路径提示。

**不**强制「先做完 tier1」；若早期 Agent 只修简单烷烃，仅通过分数与失败分布自然体现。

### 8.2 架构护栏（防错误架构锁死）

1. 静态 lint（§5）  
2. 日志强制 IUPAC 规则 ID + 涉及 Layer  
3. 连续同分且同主簇 ≥ 2 轮：prompt 注入「升层提示」——优先 L1/L2/L4 结构性修复，禁止仅在 L5 打补丁  
4. S3 边界：主规则 1 + 同簇同层附属 ≤ N-1（默认 N=3）

### 8.3 Git 门禁

- 每轮开始记录 `BASE_SHA`  
- 成功：`git add` 允许路径 → `commit`（信息含规则 ID、bench 变化）  
- 失败：`git reset --hard BASE_SHA` + `git clean` 仅限允许污路径策略（谨慎）  
- 仓库需先 `git init`（若尚未）

---

## 9. TDD 与 chem-tdd-skill

### 9.1 位置

`skills/chem-tdd-skill/`（P1，版本管理）  
编排器每轮将 skill 路径与摘要写入 cycle prompt。

### 9.2 流程（Agent 必须）

1. 从失败簇选 3–10 个 SMILES + ≥1 近邻负例  
2. 写/改 `tests/unit/...`，标记 `IUPAC: P-x`、`Layer: Lx`  
3. 跑 pytest → **必须先红**  
4. 实现 layer 代码 → 绿  
5. structure_lint → full benchmark  
6. 遵守 S3 与缓存上限  

### 9.3 与全量 bench 关系

| 门禁 | 作用 |
|------|------|
| TDD | 本轮声称的规则被精确定义且实现 |
| lint | 架构不腐化 |
| bench | 全局不回退；发现未覆盖回归 |

---

## 10. pi 调用方式

- **首选**: 非交互 CLI（print/JSON 模式）或 RPC，cwd = 项目根  
- 输入：生成的 `agent_loop/prompts/.cycle_current.md`（含本轮簇、skill、约束、命令）  
- 超时：可配（如 30–60 min/轮）  
- 失败（崩溃/超时）：本轮计失败，reset，no_improve++  

`pi_runner` 抽象接口，便于以后换 agent CLI。

---

## 11. 进度记忆格式

`agent_loop/memory/progress.md`：

```markdown
## 日志
[#hash][IUPAC P-xx.x] 标题 [+n tests, dual a%→b%]
  规则说明（≤200字）
  修改过程（100–200字）

## 进行中
（当前簇 / 阻塞）

## 其他
（≤20 行建议；被采纳后删除）
```

`STATE.json`：`iter`, `best_dual`, `last_dual`, `no_improve`, `last_cluster`, `target`, `K`。

每轮结束写 `agent_loop/memory/sessions/cycle-NNNN.json`（供 Console 侧栏）：

```json
{
  "iter": 1,
  "status": "gate_pass",
  "base_sha": "...",
  "commit_sha": "...",
  "rules": ["P-63.1.1"],
  "layers": ["L1", "L5"],
  "cluster_summary": "...",
  "prompt_path": "agent_loop/prompts/.cycle_current.md",
  "bench_before": {"dual": 0.12, "en": 0.15, "zh": 0.4, "fails": 900},
  "bench_after": {"dual": 0.14, "en": 0.17, "zh": 0.42, "fails": 880},
  "tdd": {"added": 4, "passed": true},
  "lint_ok": true,
  "decision": "commit",
  "log_excerpt": "..."
}
```

---

## 12. ChemAgent Console（网页控制台）

设计系统：`design-system/chemagent-console/`（ui-ux-pro-max）  
- 风格：**Dark Mode OLED**，密集 dashboard  
- 字体：Inter  
- 主色背景 `#0F172A`，强调/运行 `#22C55E`，危险 `#EF4444`  
- 动效低（3/10），密度高（8/10）

### 12.1 功能

| 能力 | 说明 |
|------|------|
| **Session 侧栏** | 左栏列出每轮 cycle；点选查看该轮完整上下文（簇、规则、prompt、bench、决策） |
| **启停控制** | Start / Pause（本轮结束后停）/ Stop（STOP 文件 + 可选强制 abort，需确认） |
| **Namer 试跑** | 输入 SMILES → 调本地 `namepredict` → 显示 en/zh/meta/耗时 |
| **实时日志** | 编排器与 pi 输出流式展示 |
| **Bench 视图** | dual/en/zh 趋势与失败样本表 |

### 12.2 布局（摘要）

```
Top: Logo | status | dual% | iter | Start Pause Stop
Left: Session list + fail-cluster chips
Main tabs: Overview | Cycle | Prompt | Bench | Namer | Logs
```

详见 `design-system/chemagent-console/pages/console.md`。

### 12.3 API（FastAPI）

| 方法 | 路径 | 作用 |
|------|------|------|
| GET | `/api/v1/loop/state` | 当前状态机、分数、no_improve |
| POST | `/api/v1/loop/start` | 启动/恢复循环 |
| POST | `/api/v1/loop/pause` | 本轮后暂停 |
| POST | `/api/v1/loop/stop` | 停止（body: `{force?: bool}`） |
| GET | `/api/v1/sessions` | 历史 cycle 列表 |
| GET | `/api/v1/sessions/{iter}` | 单轮上下文快照 |
| POST | `/api/v1/name` | `{smiles}` → 命名结果 |
| GET | `/api/v1/events` | SSE：`loop.state` / `log.line` / `cycle.updated` / `bench.done` |

### 12.4 与 loop 集成

- `agent_loop/events.py` 发布事件；server 订阅并推送给浏览器。  
- Start 在后台线程/任务跑 `loop.run()`，与 API 同进程共享 STATE。  
- **CLI 仍可独立跑 loop**（无 UI）；Console 为可选控制面。  
- 鉴权第一期：默认绑定 `127.0.0.1`（本机开发）；不暴露公网。

### 12.5 前端技术默认

- `web/` 静态 SPA：HTML + Tailwind（或 CDN）+ 轻量 JS；复杂度升高可迁 Vite+React。  
- 图标：Lucide SVG（禁止 emoji 图标）。  
- 图表：dual 折线 + 失败 feature 横条（少依赖，可 Chart.js）。

---

## 13. 第一期实现范围（确认后编码）

### Phase A — 基础设施

1. `git init`（若需要）、`pyproject.toml` / 依赖  
2. `tools/merge_datasets.py` → `data/merged_benchmark.json`  
3. `src/namepredict` Layer0–5 骨架 + `namer.name()`  
   - 最小规则：甲烷/乙烷/乙醇等极简可跑（证明管道）  
4. `benchmarks/benchmark.py`（分字段计分 + 归一化 + 报告）  
5. `tools/structure_lint.py`  
6. 示例 unit test + pytest 可跑  

### Phase B — 工作流

1. `agent_loop/*` 状态机、git_gate、pi_runner、events  
2. `prompts/cycle.md` 模板  
3. `skills/chem-tdd-skill/**`  
4. 单轮 dry-run（可 mock pi）与真实 pi 开关  

### Phase C — Console

1. FastAPI：`loop/sessions/name/events`  
2. `web/` Console UI（按 design-system）  
3. 启停联调 + Namer 试跑 + session 侧栏读 `cycle-*.json`  

### Phase D — 首轮验证

1. 无 pi：人工按 skill 改一轮，验证门禁  
2. 有 pi：跑 ≥1 轮自动循环，确认 commit/revert/记忆  
3. 浏览器完成：看 session → 试命名 → Pause/Stop  

**非第一期**: 完整蓝皮书规则、chebi 中文补全、性能 10ms 优化、Console 多用户鉴权。

---

## 14. 风险与缓解

| 风险 | 缓解 |
|------|------|
| 失败驱动导致局部补丁架构 | lint + 升层提示 + 禁止单 SMILES 特判 |
| chebi 过难早期 0 分 | 分字段计分；报告分 source；提升看 fails↓ |
| pi 不可用/Windows 路径 | pi_runner 可 mock；路径 pathlib |
| 函数 ≤10 行过严 | lint 强制拆分；模板展示写法 |
| 金标噪声 | 不自动改金标；可人工 quarantine 列表（后置） |
| 上下文膨胀 | 每轮新 session 或压缩；只注入簇摘要与规则指针 |
| Console 与 loop 竞态 | 单写者状态机 + 事件总线；Stop 用文件标志 |
| 误点 Stop 丢工作 | Pause 默认；Stop force 二次确认 |

---

## 15. 验收清单（设计落地后）

- [ ] 合并数据集生成且 `eval_zh` 统计合理  
- [ ] `python -m benchmarks.benchmark` 输出 en/zh/dual  
- [ ] `namer.name("CCO")` 管道通（哪怕早期不完美）  
- [ ] structure_lint 能抓超行数/错层  
- [ ] chem-tdd-skill 存在且 cycle prompt 引用  
- [ ] loop 一轮：fail→reset 与 pass→commit 两条路径可测  
- [ ] 停机条件可配置并生效  
- [ ] Console：侧栏 session、Start/Pause/Stop、Namer、SSE 日志可用  
- [ ] UI 符合 dark OLED 设计令牌与 a11y 基本项  

---

## 16. 开放小项（实现时可默认，无需再阻塞）

- 跨层 import 白名单细节  
- pi 具体 CLI 参数以本机 `pi --help` 为准  
- quarantine 金标机制延后  
- 中文名来源仅使用已有 gold，不自动翻译 chebi  
- Console 前端：静态 SPA vs Vite+React（默认静态，复杂度再升）  
- SSE vs WebSocket（默认 SSE，实现更简单）  

---

**请确认本设计**。确认后进入 `docs/superpowers/plans/` 实现计划并开始 Phase A 编码。
