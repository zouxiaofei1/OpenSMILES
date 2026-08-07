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
| 更新 Wiki | `/maintain:wiki` | 按 repowiki 流程同步文档（update / rewrite） |

> **多选执行顺序**：推荐 `dead-code → merge-files → wiki`。前两项改代码/文件结构，wiki 应最后跑以一次性反映所有改动；同时每完成一项让自动 commit 进程落盘，避免三项改动混在一个 diff 里。

> **前置约束（所有任务都遵守）**
> - **自动 commit 进程存在**：动工前先 `git status` 确认工作区状态；每轮验证通过后让自动 commit 进程（或手动 commit）落盘，避免大 diff 混在一起。
> - **测试基线非零**：当前开发中 pytest 基线有 ~576 失败（预存失败）。判回归用**相对变化**，不是绝对通过。
> - **只读文件**：绝不修改 `benchmarks/` 与 `data/merged_benchmark.json`（prompt.txt 约束）。
> - 命令不放 >15s 的 sleep。

---

## `/maintain:dead-code` — 删除死代码

### 判定标准（核心）

> **候选删除后，`benchmark` 输出不变 且 `pytest` 输出不变 → 真死代码，保留删除。**
> 任一指标变化 → 回滚（`git checkout` 该文件），不是死代码。

**基准对比指标**（前后逐项比对，全等才判死）：
- benchmark：EN acc、ZH acc、**dual acc**、**fails 数量**（report 的 `acc_dual`/`ok_*`/`fails`）
- pytest：**passed 数量**（`pytest -n auto` 汇总行）

```bash
# 基线（全量约 4070 分子，很快，一般10s内）
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time --timeout 1 2>&1 | tail -25
pytest -n auto 2>&1 | tail -5
```

> dual acc / fails 数用**取整对比**（`fails` 数必须精确相等；acc 允许最后一位浮点抖动时手工复核，`fails` 是决定性指标）。

### 候选发现方法

1. **AST 引用图 + 闭包判死**（首选）：写临时脚本在 `tools/` 下用 `ast` 解析 `src/namepredict/**`，构建"定义→引用"图，找出 `def _x` 定义了但全库无引用的函数/类。对 `_` 前缀函数尤其敏感。**引用要算闭包**——只被死函数引用的函数也是死的（沿引用图迭代直到稳定）。
2. **动态采样未覆盖**：`python tools/callchain_cli.py --n 300`（或 `--pct 100`）找"静态可达但采样从未执行"的分支——它们是死分支候选。适用"函数活着但代码路径死"的情况。
3. **已知盲区（不可只靠引用图，否则误删）**：
   - `@_register` 注册的函数：被 `ring_core_fns` 等注册表按名字动态接线，引用图不显示。删除前先查注册表消费方。
   - **转口 import**：`from namepredict.tools import side_alkyl` 这类模块级 re-export——即使模块内函数无引用，模块本身被 import 就"活着"。
   - **契约/probe 元组**：`test_l2_l3_side_facts_contract.py` 的 `_HELPERS` 等显式列出函数名的元组——删除会直接 fail 该测试，**这些函数不可删**（memory：4 函数受保护）。
   - **kind_registry 未接线 kind**：注册了但从未产生的 kind（98/157）**不一定是死**——L5 可能经其他路径派生（反例存在），删除前必须 benchmark 验证。

### 执行流程（循环多轮）

```
1. git status 确认干净起点
2. 跑基线 benchmark + pytest，记录 acc/fails/passed
3. 生成候选列表（AST 引用图 + 采样未覆盖，排除盲区）
4. 对每个候选：
   a. 删除候选（函数/import/文件）
   b. 快速定向 pytest 相关测试（如 test_l2_l3_side_facts_contract）抓明显回归
   c. 全量 pytest -n auto 对比 passed 数
   d. 全量 benchmark 对比 acc/fails
   e. 全等 → 保留删除；变化 → git checkout 回滚，标记"非死"
5. 一轮结束若删除了任何代码 → 回到步骤 2 再来一轮（删除可能暴露新的死代码）
6. 直到一轮零删除 → 完成
```

**可多次运行**：删除往往"连锁暴露"新的死代码（A 死了 → 只有 A 引用的 B 也死了），务必循环到一轮无删除才停。

### 误判陷阱（务必遵守）

- **删 layer2 `.py` 会让 contract 测试动态 glob 收缩**：`test_l2_l3_side_facts_contract` 等用 `Path.rglob("*.py")` 动态收集参数，删文件 → 参数变少 → **passed 数下降但非回归**。此时看**具体失败原因**：如果是"glob 收集的文件变少导致该测试自身断言范围变窄"，不是回归；如果出现新的失败/异常，才是回归。
- `test_layer3_functions_are_explicit`：layer3 全部函数禁止 `*args/**kwargs`——新增/改动函数注意别引入 vararg。
- 删除文件后 `grep -rln "被删模块名" src/ tests/` 必须为空（除刻意保留的契约字符串），否则有悬挂 import。

---

## `/maintain:merge-files` — 合并每层内小文件

### 触发条件

同一层（`layer0`–`layer5`、`tools/`、`layer2/scaffold/`）内存在：
- **体积 <5KB（约 <100 行）** 的文件，且
- **职责相似**（同属一个功能簇：如同层多个 `side_*` 拓扑探针、多个 `*_names` 命名表、多个 scaffold 薄层）。

### 流程

```
1. git status 确认
2. 列层内文件 + 字节数/行数：ls -la + wc -l
3. 按职责聚类小文件，形成合并方案（例：layer3 的 side_facts/side_alkyl/side_alkoxy 属"拓扑事实簇"；alkyl_sys_names/alkoxy_names/amino_side 属"命名表簇"）
4. 用 AskUserQuestion 让用户确认每个合并组（可多选跳过）
5. 执行合并：
   a. 把目标文件内容并入宿主文件（宿主通常是层内职责最广者）
   b. 全库 grep -rln "被并入模块名" src/ tests/ 更新所有 import 路径
   c. 删除被并入的源文件
   d. 处理 __all__ / __init__ re-export（如有）
6. 验证零回归：
   - pytest -n auto 对比 passed（注意 contract glob 收缩陷阱，见上）
   - python -m benchmarks.benchmark_parallel ... 对比 acc/fails
   - test_l2_l3_side_facts_contract 的 LAYER3.rglob 扫描：合并后 make_match 等显式参数契约必须仍满足
7. 提交
```

### 约束

- 合并后**单文件 ≤1000 行**（prompt.txt 硬约束）。
- **跨层不合并**——`tools/` 是层无关工具，只与 `tools/` 内合并；`layer3/leaves/` 与 `layer3/` 其余文件不混。
- **契约测试保护的文件慎重**：`leaves/protocol.py`（make_match 显式参数契约）、`side_facts.py`（test_l2_l3_side_facts_contract 按路径扫描）——合并它们要同步更新测试的 `LAYER3`/`LEAF_PROTOCOL` 路径常量，否则测试立即 fail。
- 合并会减少文件数 → contract 测试 glob 参数收缩 → passed 下降是预期的，看具体失败而非数字。

---

## `/maintain:wiki` — 更新 RepoWiki

复用 `.claude/skills/repowiki/SKILL.md` 的流程，按重构规模选粒度：

| 场景 | 用 |
|------|-----|
| 小改动（改了几个函数/常量） | **update**：`git diff --name-only HEAD` → 与 `docs/wiki/log.md` 比对 → 只重写受影响页面 |
| 大重构（文件搬移、模块删除、机制替换） | **rewrite 受影响页面**：重写 layer3-substituents.md 等架构页 + 修正其他页面对已删模块的引用 |
| 链接/孤立检查 | `lint`：孤立页面、失效源码锚点、过时页面 |

### 标准流程

```
1. git diff --name-only HEAD  确定变更文件
2. 判定重构规模 → update or rewrite
3. 扫描受影响页面，修正：
   - 文件清单（新增/删除/搬移的模块，行数）
   - 源码锚点 file:line（删行后行号偏移，逐页核对）
   - 已删除模块的残留引用（grep -rn "已删模块名" docs/wiki）
   - mermaid 数据流图（结构变化必须同步）
4. 更新 index.md 统计（src .py 总数 / 总行数：find src -name "*.py" | wc -l）
5. 追加 log.md 记录（append-only 表格）
6. 验证：mermaid ``` 配对、无已删模块引用残留
```

### 经验（2026-08 重构示例）

- 拓扑事实层 `side_facts/side_alkyl/side_alkoxy/aryl_sub/aryl_depth2/leaves` 从 layer2/tools 迁入 layer3 后，**layer2 页面里"L2 负责侧链分类"的说法反转**——L2 只反向借用少量函数。搜 layer2 页面里所有 `layer2/side_*`、`layer2/aryl_*`、`layer2/leaves` 引用逐一改路径。
- 删除模块（如 `ring_namer`/`heteroaryl_sub`/`cycloalkyl_names`）后，全 wiki `grep -rn` 这些名字，除了"已删除"说明外不应有残留；guides/overview 的目录树里也常藏着。
- 重构后 `layer3-substituents.md` 这类"机制中心"页面优先全量 rewrite，而不是逐段 edit——它的提取算法、数据流图、文件清单全变了。

---

## 附：常用命令速查

```bash
# 测试（多线程）
pytest -n auto
# 定向契约测试
pytest tests/unit/test_l2_l3_side_facts_contract.py -v
# benchmark 全量
python -m benchmarks.benchmark_parallel --data data/merged_benchmark.json --time --timeout 1
# 动态调用链（找采样未覆盖）
python tools/callchain_cli.py --n 300
# 静态调用图（layer2）
python tools/call_graph.py
# 文件体积
find src/namepredict -name "*.py" | xargs wc -l | sort -n | head -30
# 更新 index 统计
echo -n "src .py 数: "; find src -name "*.py" -not -path "*/__pycache__/*" | wc -l
echo -n "src 总行数: "; find src -name "*.py" -not -path "*/__pycache__/*" -exec cat {} + | wc -l
```
