---
name: chem-tdd-skill
description: SMILES→IUPAC 编写的 TDD 流程
---

# chem-tdd-skill

本 skill 约束 **每一轮** Agent 对 `namepredict` 的规则实现。编排器会将本路径注入 cycle prompt；**必须全文遵守**，不得跳过红→绿门禁。

## 目标

在不破坏架构与全局准确率的前提下，用 **可复现的单元测试** 精确定义本轮 IUPAC 规则，再实现 Layer 代码，使该规则通过 TDD，并通过 lint 与全量 benchmark 门禁。

## 硬流程（先测后码；先红后绿）

**禁止**先写/改生产代码再补测试。每轮必须按序：

1. **读上下文**：本轮失败簇摘要、建议 Layer、相关 `docs/cleaned` / `docs/iupac` 路径、`agent_loop/memory/progress.md`。
2. **抽例**：从失败簇抽取 **3–10 个 SMILES** 作为正例，并至少 **1 个近邻负例**（结构相近但应得到不同名称，或本规则不应误伤的已正确分子）。
3. **写测试**（仅 `tests/unit/`）：按模板新建/扩展用例；文件头标注 `IUPAC:` 与 `Layer:`。
4. **跑 pytest → 必须先红**：在尚未实现本轮规则时，新增断言必须失败。若已绿，说明用例无效或规则已存在——重写用例，禁止「假绿」。
5. **实现代码**（仅允许路径内 Layer/入口/cache）：使上述用例变绿；遵守 S3 与行数/函数长度约束。
6. **再跑 pytest → 必须绿**。
7. **structure_lint** → **full benchmark**。
8. 按门禁决定 commit 或依赖编排器 reset；写 **回合小结**。

红→绿证据：同一批新增/修改的测试，实现前失败、实现后通过。不得删除/弱化断言来换绿。

## 抽例规则（3–10 SMILES + ≥1 负例）

| 类型 | 数量 | 要求 |
|------|------|------|
| 正例 | 3–10 | 来自本轮失败簇；覆盖规则边界（最短链/最长相关、有无侧链、关键异构等） |
| 负例 | ≥1 | 近邻但 **不应** 被本规则改名错误；或应保持现有正确输出 |

- 负例的期望名取自金标或当前已正确行为；断言其输出 **不被** 本轮改动破坏。
- 禁止只用 1–2 个正例「过关」；少于 3 正例或 0 负例视为流程违规。
- 用例数据写在测试文件的 `CASES` 中，不写进生产代码。

## 双语断言规则

- 调用：`SMILESNNamer().name(smiles)`，先 `assert r.success`。
- **英文**：`assert normalize_en(r.en) == normalize_en(en)`（轻度归一化：strip、小写、空白合并）。
- **中文**：当 `zh is not None` 时  
  `assert normalize_zh(r.zh) == normalize_zh(zh)`（strip 后严格相等；不转繁简、不忽略标点）。
- 金标无中文时：`zh=None`，**跳过**中文断言，但仍须断言英文（若本簇测英）。
- 不得手写「近似相等」、子串包含、忽略位次等放宽匹配。
- 导入：`from namepredict.constants import normalize_en, normalize_zh`。

## 测试文件头（强制）

每个本轮新增/主改的 `tests/unit/test_*.py` 顶部注释必须含：

```text
# IUPAC: P-XX.X.X
# Layer: Lx
```

- `IUPAC:` 为本轮主规则 ID（蓝皮书条款，如 `P-63.1.1`）；附属规则在小结中列出。
- `Layer:` 为主要改动层，如 `L2` 或 `L1,L5`（逗号分隔，无空格或单空格均可，但须可 grep）。
- 可参考模板：`skills/chem-tdd-skill/templates/test_rule.py.tmpl`。

## S3 边界（N=3）

- **默认**：每轮 **1 个主规则**。
- **上限**：主规则 1 + 同簇、同层附属修复 ≤ **N−1**，默认 **N=3**，即同簇同层合计改动规则点 **≤ 3**。
- 超出范围的失败留给后续轮次；禁止「顺手修一片」。
- 连续同分且同主簇 ≥ 2 轮：优先 L1/L2/L4 结构性修复，**禁止**仅在 L5 打补丁。

## 禁止事项

1. **禁止** `if smiles == "..."` / 哈希表按完整 SMILES 映射名称 等 **单 SMILES 特判**。
2. **禁止** 修改 benchmark **金标**（`data/*`、`chebi20_test_1k.json`、`smiles_tiers.json`、`data/merged_benchmark.json`）及 `benchmarks/benchmark.py` 计分逻辑。
3. **禁止** 用深度学习/LLM 做命名；俗名+例外 **cache ≤ 100**。
4. **禁止** 居所取代型投机捷径；须按 IUPAC 系统逻辑实现。
5. **禁止** 先改生产代码再补测、或改测迁就错误实现。
6. **禁止** 为过门禁而删除/跳过失败测试、放宽归一化、改 lint 白名单（除非人类任务明确要求）。
7. 生产代码仅可写：`src/namepredict/` 入口、`layer0`–`layer5`、受限 `cache`；测试写 `tests/unit/`。  
   单文件 ≤ 500 行；单函数/方法体 ≤ 10 行。

## 命令（逐条执行）

在项目根目录：

```bash
# TDD：先红后绿
pytest tests/unit -q

# 架构 lint
python tools/structure_lint.py

# 全量（或编排器指定的）benchmark
python -m benchmarks.benchmark --data data/merged_benchmark.json
# 调试可加 --limit N；门禁轮次以编排器 {{bench_cmd}} 为准
```

建议顺序：pytest（红）→ 实现 → pytest（绿）→ structure_lint → benchmark。

## 可写 / 只读路径

| 可写 | 只读 |
|------|------|
| `src/namepredict/`（namer 入口、layer0–5、cache） | `data/*` 金标与合并集 |
| `tests/unit/` | `benchmarks/benchmark.py` |
| （记忆区由编排器为主）`agent_loop/memory/` | `docs/iupac/**`、`docs/cleaned/**` |

## 门禁含义（Agent 自检）

| 门禁 | 作用 |
|------|------|
| TDD (pytest) | 本轮规则被测试精确定义且实现为真 |
| structure_lint | 层职责、行数、无 SMILES 特判、缓存上限 |
| benchmark | 全局 dual 不回退 >0.5%；dual↑ 或 fails↓ 才算提升 |

## 回合小结格式（每轮结束必须输出）

```markdown
## 回合小结
- **iter**: <N>
- **IUPAC**: <主规则ID>；（附属: ...）
- **Layer**: <Lx,...>
- **簇**: <cluster key / 一句话>
- **TDD**: 正例 <n> + 负例 <m>；红→绿：是/否；新增测试文件：<paths>
- **改动文件**: <list>
- **pytest**: pass/fail
- **lint**: pass/fail
- **bench**: dual <before>% → <after>%；fails <a> → <b>
- **决策建议**: commit / reset（及原因）
- **规则说明**（≤200字）: ...
- **修改过程**（100–200字）: ...
```

进度记忆（若本轮负责写 `progress.md` 日志行）对齐：

```text
[#hash][IUPAC P-xx.x] 标题 [+n tests, dual a%→b%]
  规则说明（≤200字）
  修改过程（100–200字）
```

## 检查清单

开始实现前与提交前，对照 `skills/chem-tdd-skill/checklist.md` 逐项打勾。

## 快速参考

- 模板：`skills/chem-tdd-skill/templates/test_rule.py.tmpl`
- 失败簇工具：`tools/fail_cluster.py` → `cluster_failures`
- 归一化：`namepredict.constants.normalize_en` / `normalize_zh`
- 入口：`namepredict.namer.SMILESNNamer`
