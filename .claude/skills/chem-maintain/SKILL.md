---
name: chem-maintain
description: namepredict 项目维护三件套——(1) 删除死代码：以 benchmark + pytest 输出不变判死，可循环多轮；(2) 合并每层内 <5KB、职责相似的过小文件；(3) 更新 RepoWiki 反映重构。触发词：维护、死代码、删死代码、清理、合并文件、小文件合并、更新wiki、维护任务、maintain、dead code。
---

# chem-maintain — 项目维护三件套

为 `namepredict`（SMILES → IUPAC 双语命名引擎）执行三类维护任务。**skill 被调用后，第一件事用 `AskUserQuestion`（multiSelect=true）询问用户执行哪一项**，按选择执行：

| 选项 | 子命令 | 说明 |
|------|--------|------|
| 删除死代码 | `/maintain:dead-code` | AST/采样找候选 → 删 → benchmark+pytest 输出不变=真死 → 循环多轮 |
| 合并小文件 | `/maintain:merge-files` | 每层内 <5KB 且职责相似的文件合并，更新 import |
| 更新 Wiki | `/maintain:wiki` | 按 repowiki 流程同步文档（update / rewrite / generate） |

> **多选执行顺序**：推荐 `dead-code → merge-files → wiki`。前两项改代码/文件结构，wiki 应最后跑以一次性反映所有改动；同时每完成一项就手动 `git add` + commit 落盘，避免三项改动混在一个 diff 里。

> **前置约束（所有任务都遵守）**
> - **手动 commit**：当前仓库无自动 commit hook（.claude/settings 仅放行 `git add *` / `git commit` 权限）。动工前先 `git status` 确认工作区；每轮验证通过后手动 commit 落盘，避免大 diff 混在一起。
> - **测试基线全绿（2026-09 实测）**：`pytest -n auto` = **1702 passed / 0 failed**。判回归用**绝对标准**——出现任何失败即回归；无需相对 passed 数对比。定向测试见"契约测试雷达"。
> - **benchmark 基线**：`data/merged_benchmark.json` = **3990 分子，实测 ~14s**（19 worker）。summary 行 `en=.. zh=.. dual=.. fails=N`。字段 `acc_en/acc_zh/acc_dual/ok_*/fails` 均在。**fails 数精确相等是决定性指标**，acc 允许最后一位浮点抖动时手工复核。
> - **只读文件**：绝不修改 `benchmarks/` 与 `data/merged_benchmark.json`（prompt.txt 约束）。
> - 命令不放 >15s 的 sleep。
> - **根 `tools/` ≠ 源码**：根 `tools/`、`tmp/` 整个被 gitignore（prompt.txt"临时文件写 tools"），是分析/临时脚本区，非版本化源码。**源码工具在 `src/namepredict/tools/`**（anchored_table / chain / block_cut / free_to_yl，被 git 跟踪）。别混淆两处。

---

## `/maintain:dead-code` — 删除死代码

### 判定标准（核心）

> **候选删除后，`benchmark` 输出不变 且 `pytest` 不变（0 失败） → 真死代码，保留删除。**
> 任一指标变化 → 回滚（`git checkout` 该文件），不是死代码。

**基准对比指标**（前后逐项比对，全等才判死）：
- benchmark：EN acc、ZH acc、dual acc、**fails 数量**（summary 行的 `en=.. zh=.. dual=.. fails=N`；`fails` 精确相等才判过）
- pytest：**0 failed**（`pytest -n auto` 汇总行；删除会触发"参数化 glob 收缩"导致的 passed 下降，见误判陷阱——不视为回归）

```bash
# 基线（全量 3990 分子，实测 ~14s）
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time --timeout 1 2>&1 | tail -6
pytest -n auto 2>&1 | tail -5
```

### 候选发现方法

1. **AST 引用图 + 闭包判死**（首选）：写临时脚本在根 `tools/` 下用 `ast` 解析 `src/namepredict/**`，构建"定义→引用"图，找出 `def _x` 定义了但全库无引用的函数/类。对 `_` 前缀函数尤其敏感。**引用要算闭包**——只被死函数引用的函数也是死的（沿引用图迭代直到稳定）。
2. **动态采样未覆盖**：`python tools/callchain_cli.py --n 300`（或 `--pct 100`）找"静态可达但采样从未执行"的分支。**注意此脚本在根 `tools/`（gitignore 临时区，非版本化）**——文件还在时可用；若已被清理则退化为纯 AST 法。
3. **已知盲区（不可只靠引用图，否则误删）**：
   - **数据表注册型（不是函数引用）**：kind/scaffold/词干通过注册表按 id 接线，AST 引用图不显示。典型：`layer1/fg_registry.FG_SPECS`（FgSpec 表，principal / functional_group_inventory 等派生唯一来源）、`layer2/ring_scaffold` 的 ScaffoldSpec/registry（kind_registry 启动时经 `_load_from_scaffold_specs` 从 `all_specs` 回填 KindMeta）——删注册表条目会连带词干/L5，须 pytest（`test_scaffold_registry_single_source` / `test_ring_scaffold_expression_contract` / `test_retained_registry`）+ benchmark 验证。
   - **转口 import**：模块级 re-export——即使模块内函数无引用，模块本身被 `from namepredict.layerX import ...` import 就"活着"。删前查各层 `__init__.py` 与 namer.py 的 re-export。
   - **契约名单 / AST 扫描（挂名的"4 函数"已过时）**：`tests/unit/test_l2_parent_core_contract.py:22` 的 `_HELPERS = {_parent_dict, _longest_from, _longest_chain}`（**3 个，非 4**）断言 L2 producer 不得从 `parent_selector` import 它们——三个 helper 现已迁至 `tools/chain.py`、`layer2/chain_walk.py`、`layer2/principal_expression.py`，**改动/重排这几个 helper 所在模块要小心此契约**。`test_l2_l3_side_facts_contract.py` 名字里的 side_facts/leaves 已名不副实，现用 AST 扫描 `tools/chain.py` 的拓扑约束（无 dict 分发、无命名依赖、dataclass 只暴露拓扑字段）——**别被旧测试名误导**。
   - **kind 注册但未产生，不一定是死**：kind_registry 现在是 KindMeta 元数据壳，词干全来自 ring_scaffold 回填；旧 `ring_producers`/`unsat_producers` 引导链已删（见 layer2/candidates.py docstring）。任何"kind 存在但命名路径没走到"的反例都要以 benchmark 为准，删除前必须验证。

### 执行流程（循环多轮）

```
1. git status 确认干净起点（有未提交改动先 commit）
2. 跑基线 benchmark + pytest，记录 acc/fails / 0-failed
3. 生成候选列表（AST 引用图 + 采样未覆盖，排除盲区）
4. 对每个候选：
   a. 删除候选（函数/import/文件）
   b. 快速定向 pytest 契约雷达里的相关测试抓明显回归
   c. 全量 pytest -n auto（0 failed 即过；passed 因 glob 收缩下降属预期）
   d. 全量 benchmark 对比 acc/fails（fails 精确相等）
   e. 全等 → 保留删除；变化 → git checkout 回滚，标记"非死"
5. 一轮结束若删除了任何代码 → 回到步骤 2 再来一轮（删除可能暴露新的死代码）
6. 直到一轮零删除 → 完成
```

**可多次运行**：删除往往"连锁暴露"新的死代码（A 死了 → 只有 A 引用的 B 也死了），务必循环到一轮无删除才停。

### 误判陷阱（务必遵守）

- **删 layer2 `.py` 会让 `test_l2_parent_core_contract` 的 rglob 参数化收缩**：`_L2.rglob("*.py")` 在 import 时求值喂给 `@pytest.mark.parametrize`，删文件 → 参数变少 → **passed 数下降但非回归**。此时看**具体失败原因**：出现新的失败/异常才是回归；只是该测试断言范围变窄不是。`test_l2_l3_side_facts_contract` **无参数化**（函数体内 rglob 扫描），删 layer3 文件不降 passed、只缩扫描范围，同理看失败而非数字。
- **layer3 全函数禁止 `*args/**kwargs`**：约束在 `test_l2_l3_side_facts_contract.py::test_layer3_functions_are_explicit`（rglob 扫 layer3 全函数；不是独立测试文件）——新增/改动 layer3 函数别引入 vararg。
- **`tools/chain.py` 被一长串 AST 契约钉死**：`test_l2_l3_side_facts_contract.py` 里 `test_side_facts_*` 断言它公开函数带类型、无 dict 分发/命名 token 表/name_leaf 依赖、dataclass 仅拓扑字段——改动 chain.py 会让一串测试 fail（测试名旧，测试本身活着）。
- 删除文件后 `grep -rln "被删模块名" src/ tests/` 必须为空（除契约字符串、`docs/wiki/log.md` 历史），否则有悬挂 import。

### 契约测试雷达（dead-code / merge-files 的定向回归网）

| 测试 | 保护对象 | 删除/改动触发 |
|------|---------|--------------|
| `test_l2_parent_core_contract.py` | `_HELPERS`3 个 helper 不得从 parent_selector import；`_L2.rglob` 参数化 | 删 L2 .py → 参数化收缩 |
| `test_l2_l3_side_facts_contract.py` | `tools/chain.py` 拓扑约束 + layer3 无 vararg | 改 chain.py → 一串 fail；新 layer3 vararg → fail |
| `test_scaffold_registry_single_source.py` | ring_scaffold 词干/retained 是 kind_registry 唯一来源 | 删/改 scaffold 注册条目 |
| `test_ring_scaffold_expression_contract.py` / `test_retained_registry.py` | scaffold_identity 解析、retained 匹配 | 动 ring_scaffold 匹配逻辑 |
| `test_l5_no_layer2_private.py` | layer5 不得 import layer2（`_L5.glob`） | 给 layer5 加 layer2 import |

---

## `/maintain:merge-files` — 合并每层内小文件

> **2026-09 现状**：layer 已扁平化——`layer2/scaffold/`、`layer3/leaves/` 等子目录全部删除；`side_facts.py`/`side_alkyl.py`/`side_alkoxy.py`/`aryl_*`/`alkyl_sys_names`/`alkoxy_names`/`leaves/protocol.py` 等旧小文件均已并入/移除。当前各层文件普遍 >40 行，**此子命令可选对象很少**。仅在新文件堆积出 <5KB 同职责簇时才用，勿为合并而合并。

### 触发条件

同一层（`layer0`–`layer5`、`src/namepredict/tools/`）内存在：
- **体积 <5KB（约 <100 行）** 的文件，且
- **职责相似**（同属一个功能簇：如 layer3 的 amino_side/as_substituent 薄封装、tools 下的小原语）。

### 流程

```
1. git status 确认（未提交改动先 commit）
2. 列层内文件 + 字节数/行数：ls -la + wc -l
3. 按职责聚类小文件，形成合并方案
4. 用 AskUserQuestion 让用户确认每个合并组（可多选跳过）
5. 执行合并：
   a. 把目标文件内容并入宿主文件（宿主通常是层内职责最广者）
   b. 全库 grep -rln "被并入模块名" src/ tests/ 更新所有 import 路径
   c. 删除被并入的源文件
   d. 处理 __all__ / __init__ re-export（如有）
6. 验证零回归（跑上表"契约测试雷达"全部 + benchmark）：
   - pytest -n auto：0 failed（删 L2 文件会让 test_l2_parent_core_contract 参数化收缩，看失败而非数字）
   - python -m benchmarks.benchmark_parallel ... 对比 acc/fails（fails 精确相等）
7. 提交
```

### 约束

- 合并后**单文件 ≤1000 行**（prompt.txt 硬约束）。
- **跨层不合并**——`src/namepredict/tools/` 是层无关工具，只与其内部合并；不与 layer3 混。
- **契约测试钉死的文件慎重**：`src/namepredict/tools/chain.py`（`test_l2_l3_side_facts_contract` 的 AST 约束）、`layer2/ring_scaffold.py`（词干唯一来源，`test_scaffold_registry_single_source`）——合并它们会立即触发雷达测试 fail。
- 合并减少文件数 → `test_l2_parent_core_contract` 的 rglob 参数化参数收缩 → passed 下降是预期的，看具体失败而非数字。

---

## `/maintain:wiki` — 更新 RepoWiki

复用 `.claude/skills/repowiki/SKILL.md` 的流程。**repowiki 现提供 7 个子命令**（非早期三粒度）：`generate`（全量/两阶段编排生成）、`update`、`modify <page>`、`supplement <page>`、`rewrite <page>`、`plan`、`lint`。按重构规模选粒度：

| 场景 | 用 |
|------|-----|
| 小改动（改了几个函数/常量） | **update**：`git diff --name-only HEAD` → 与 `docs/wiki/log.md`（append-only 表，截至 2026-09-05）比对 → 只重写受影响页面 |
| 新增/补页 | **generate** 或 **supplement \<page\>** |
| 大重构（文件搬移、模块删除、机制替换） | **rewrite 受影响页面**：重写 `architecture/layer3-substituents.md` 等架构页 + 修正其他页面对已删模块的引用 |
| 链接/孤立检查 | **lint**：孤立页面、失效源码锚点、过时页面 |

### 标准流程

```
1. git diff --name-only HEAD  确定变更文件
2. 判定重构规模 → update / rewrite / generate
3. 扫描受影响页面，修正：
   - 文件清单（新增/删除/搬移的模块，行数）
   - 源码锚点 file:line（删行后行号偏移，逐页核对）
   - 已删除模块的残留引用（grep -rn "已删模块名" docs/wiki —— 除 log.md 历史与"已删"脚注外不应有）
   - mermaid 数据流图（结构变化必须同步）
4. 更新 `index.md` 统计（当前口径：含 tools/ 与根目录、排除 __pycache__/cache，2026-09 为 68 .py / 9,853 行 / 15 页）
5. 追加 `log.md` 记录（append-only 表格，格式 `| 时间 | 操作 | 涉及页面 | 触发者 |`）
6. 验证：mermaid ``` 配对、无已删模块引用残留、repowiki lint
```

### 现状与注意（2026-09 核对）

- wiki 已较新（截至 2026-09-05）：`layer3-substituents.md` 等架构页存在；旧模块残留引用已基本清零——`guides/overview.md` **已不存在**（目录树实际在 `architecture/overview.md` 的"文件组织"节，且已剔除已删模块）。删模块后仍要 `grep -rn` 全 wiki 防新残留，但不要再指望 overview 目录树里藏旧引用。
- `repowiki/SKILL.md` 自身"文档结构"树里的层 3/4 文件名有误（写 `layer3-numbering.md`/`layer4-locants.md`，实际为 `layer3-substituents.md`/`layer4-numbering.md`）——若按 repowiki 流程撞上，以实际文件为准，属 repowiki skill 待修的 bug。
- log.md 的早期条目（ring_namer/heteroaryl_sub/cycloalkyl_names/side_* 迁移等）是**历史记录**，勿当残留清理。

---

## 附：常用命令速查

```bash
# 测试（多线程，全绿基线 1702 passed / 0 failed）
pytest -n auto
# 定向契约测试（dead-code/merge 回归雷达）
pytest tests/unit/test_l2_parent_core_contract.py tests/unit/test_l2_l3_side_facts_contract.py -v
pytest tests/unit/test_scaffold_registry_single_source.py tests/unit/test_l5_no_layer2_private.py -v
# benchmark 全量（3990 分子 ~14s，fails 精确相等判过）
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time --timeout 1
# 动态调用链（找采样未覆盖；脚本在 gitignore 的根 tools/，非版本化）
python tools/callchain_cli.py --n 300
# 静态调用图（layer2；同上在根 tools/）
python tools/call_graph.py
# 源码文件体积
find src/namepredict -name "*.py" -not -path "*/__pycache__/*" | xargs wc -l | sort -n | head -30
# 更新 index 统计（口径与 docs/wiki/index.md 一致：含根 tools/，排除 __pycache__/cache）
find . -name "*.py" -not -path "*/__pycache__/*" -not -path "./.venv/*" -not -path "*/cache/*" | wc -l
```
