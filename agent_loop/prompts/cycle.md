# Agent Cycle — 第 {{iter}} 轮

你是 SMILES→IUPAC 规则命名系统（`namepredict`）的实现 Agent。本轮由编排器注入上下文；**必须**按 TDD 与硬约束完成本轮改动。

## 本轮上下文

- **迭代编号**: `{{iter}}`
- **当前 dual 准确率**: `{{dual}}`
- **失败簇摘要**:

```
{{cluster}}
```

- **TDD Skill 路径**: `{{skill_path}}`（请先阅读并全文遵守）
- **额外约束**:

```
{{constraints}}
```

- **门禁 benchmark 命令**（最终以该命令为准）:

```
{{bench_cmd}}
```

## 硬约束（违反即本轮作废）

1. **禁止 ML / LLM 命名**：不得用深度学习或调用外部 LLM 生成化学名；俗名+例外 cache **≤ 100**。
2. **禁止 SMILES 特判**：禁止 `if smiles == "..."`、完整 SMILES→名称哈希表等单分子捷径。
3. **禁止改金标与计分**：不得修改 `data/*`、`chebi20_test_1k.json`、`smiles_tiers.json`、`data/merged_benchmark.json` 及 `benchmarks/benchmark.py` 计分逻辑。
4. **行数限制**：单文件 ≤ 500 行；单函数/方法体 ≤ 10 行。
5. **只改允许路径**（其它路径只读）：
   - `src/namepredict/`（入口、layer0–layer5、受限 cache）
   - `tests/unit/`
   - 必要时 `agent_loop/memory/` 会话/进度（若编排器要求写小结到指定处）
6. **Windows 兼容**：路径用 pathlib 思维；勿写死 Unix-only 假设。

## S3 粒度（N=3）

- **默认**：本轮只实现 **1 个主规则**。
- **上限**：主规则 1 + 同簇、同层附属修复 ≤ **N−1**，合计规则点 **≤ 3**。
- 超出范围的失败留给后续轮次；禁止「顺手修一片」。
- 若连续同分且同主簇 ≥ 2 轮：优先 L1/L2/L4 结构性修复，**禁止**仅在 L5 打补丁。

## 强制 TDD 流程

严格遵循 `{{skill_path}}`（`skills/chem-tdd-skill/`）：

1. 读失败簇与相关 docs；从簇中抽 **3–10 正例 + ≥1 近邻负例**。
2. **先写测试**（`tests/unit/`），文件头标注 `# IUPAC:` 与 `# Layer:`。
3. **pytest 必须先红**（未实现时断言失败）；禁止假绿。
4. **再实现** Layer/入口代码使测试变绿。
5. **pytest 再绿** → `python tools/structure_lint.py` → 运行 `{{bench_cmd}}`。
6. 不得删除/弱化断言、放宽归一化来换绿。

双语断言：`normalize_en` / `normalize_zh`；金标无中文时跳过中文断言。

## 建议命令顺序

```bash
pytest tests/unit -q          # 红
# 实现规则
pytest tests/unit -q          # 绿
python tools/structure_lint.py
{{bench_cmd}}
```

## 完成后自列小结

回合结束时，用简短列表写清：

1. **主规则**（IUPAC ID）与 **Layer**
2. **改动文件**列表（相对仓库根）
3. **新增/修改测试**与抽例数量（正例/负例）
4. **红→绿证据**（实现前失败、实现后通过）
5. **structure_lint / pytest / benchmark** 结果摘要
6. **未覆盖失败**（留给后续轮次）
7. **风险或 follow-up**

不要提交 git（由编排器 GitGate 决定 commit 或 reset）。
