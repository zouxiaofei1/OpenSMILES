---
name: chem-tdd-skill
description: SMILES→IUPAC 编写的 TDD 流程
---

# chem-tdd-skill

本 skill 包含  Agent 编写 `namepredict` 代码时 TDD的规则实现。

## 目标

在不破坏架构与全局准确率的前提下，用 **可复现的单元测试** 精确定义本轮 IUPAC 规则，再实现 Layer 代码，使该规则通过 TDD，并通过 lint 与全量 benchmark 门禁。

## TDD流程（按照本流程执行）

1. **抽例**：从失败簇抽取 **3–5 个 SMILES** 作为正例，并至少 **1 个近邻负例**（结构相近但应得到不同名称，或本规则不应误伤的已正确分子）。
2. **写测试**（仅 `tests/unit/`）：按模板新建/扩展用例；文件头标注 `IUPAC:` 与 `Layer:`。
3. **跑 pytest → 必须先红**：在尚未实现本轮规则时，新增断言必须失败。若已绿，说明用例无效或规则已存在——重写用例，禁止「假绿」。
4. **实现代码**（仅允许路径内 Layer/入口/cache）：使上述用例变绿；遵守 S3 与行数/函数长度约束。
5. **再跑 pytest → 必须绿**。


红→绿证据：同一批新增/修改的测试，实现前失败、实现后通过。不得删除/弱化断言来换绿。


## 双语断言规则

- 调用：`SMILESNNamer().name(smiles)`，先 `assert r.success`。
- **英文**：`assert normalize_en(r.en) == normalize_en(en)`（轻度归一化：strip、小写、空白合并）。
- **中文**：当 `zh is not None` 时  
  `assert normalize_zh(r.zh) == normalize_zh(zh)`（strip 后严格相等；不转繁简、不忽略标点）。
- 金标无中文时：`zh=None`，**跳过**中文断言，但仍须断言英文（若本簇测英）。
- 不得手写「近似相等」、子串包含、忽略位次等放宽匹配。
- 导入：`from namepredict.constants import normalize_en, normalize_zh`。

## 测试文件头

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

## 门禁含义（Agent 自检）

| 门禁 | 作用 |
|------|------|
| TDD (pytest) | 本轮规则被测试精确定义且实现为真 |
| structure_lint | 层职责、行数、无 SMILES 特判、缓存上限 |
| benchmark | 全局 dual 不回退 >0.5%；dual↑ 或 fails↓ 才算提升 |
