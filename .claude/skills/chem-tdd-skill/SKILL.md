---
name: chem-tdd-skill
description: SMILES→IUPAC 编写的 TDD 流程
---

# chem-tdd-skill

本 skill 包含  Agent 编写 `namepredict` 代码时 TDD的规则实现。

## 目标

TDD 和 benchmark 一起作为代码质量的双重门禁。TDD用以保障以实现功能的可靠性。

## TDD流程（按照本流程执行）

1. **抽例**：从失败簇抽取 **3–5 个 SMILES** 作为正例，并至少 **1 个近邻负例**（结构相近但应得到不同名称，或本规则不应误伤的已正确分子）。
2. **写测试**（仅 `tests/unit/`）：按模板新建/扩展用例；文件头标注 `IUPAC:` 与 `Layer:`。
3. **跑 pytest → 必须先红**：在尚未实现本轮规则时，新增断言必须失败。若已绿，说明用例无效或规则已存在——重写用例，禁止「假绿」。
4. **实现代码**（仅允许路径内 Layer/入口/cache）：使上述用例变绿；遵守行数/函数长度约束。
5. **再跑 pytest → 必须绿**。


## 测试文件头格式

每个本轮新增/主改的 `tests/unit/test_*.py` 顶部注释包含：

```text
# IUPAC: P-XX.X.X
# Layer: Lx
```

- `IUPAC:` 为本轮主规则 ID（蓝皮书条款，如 `P-63.1.1`）；附属规则在小结中列出。
- `Layer:` 为主要改动层，如 `L2` 或 `L1,L5`（逗号分隔，无空格或单空格均可，但须可 grep）。

## 禁止事项

1. **禁止** `if smiles == "..."` / 哈希表按完整 SMILES 映射名称 等 **单 SMILES 特判**。
2. **禁止** 未经明确允许修改 benchmark数据及 `benchmarks/benchmark.py` 计分逻辑。
3. **禁止** 用深度学习/LLM 做命名；俗名+例外 **cache ≤ 100**。
4. **禁止** 为过门禁而删除/跳过失败测试、放宽归一化、改 lint 白名单（除非人类任务明确要求）。


## 命令

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


## 其他
红→绿证据：同一批新增/修改的测试，实现前失败、实现后通过。不得删除/弱化断言来换绿。
生产代码仅可写：`src/namepredict/` 入口、`layer0`–`layer5`、受限 `cache`；测试写 `tests/unit/`。  