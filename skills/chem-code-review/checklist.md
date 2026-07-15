# chem-code-review 检查清单

检查 Agent 在 **开始审查时** 与 **给出裁决前** 各过一遍。任一项未满足不得 **PASS**。

## A. 选题（与代码解耦；先做）

- [ ] 已读本轮意图、`workstate.md` 候选/日志
- [ ] 已用「选题质疑清单」自问（过窄、错层、低收益、依赖未就绪、特判气质等）
- [ ] 选题同意 **或** 已准备 `CHALLENGE`（反对点 + 替代选题 + 原题条件 + diff 处置）
- [ ] 未因「代码很好」跳过选题质疑

## B. 范围

- [ ] 已取得改动文件列表（`git diff` / 实现方报告）；无 diff 则标「仅选题」
- [ ] 仅 `src/namepredict/`（入口、layer0–5、cache）与 `tests/unit/`（及约定记忆/workstate）
- [ ] 未改 `data/*`、`data/merged_benchmark.json`、`benchmarks/benchmark.py`

## C. 自动化门禁（有实现 diff 时）

- [ ] `python tools/structure_lint.py --root src/namepredict` → ok
- [ ] `pytest tests/unit -q` → 绿（或明确无测改动且说明理由）
- [ ] `python -m benchmarks.benchmark --data data/merged_benchmark.json`
- [ ] 记录 **中英文双重准确: X%**；与改前对比
- [ ] dual 回退 **未** >0.5%（否则 ROLLBACK）

## D. 体量与函数式

- [ ] 每个改动函数体非空行 ≤10；超限已要求拆分
- [ ] 单文件 ≤1000 行
- [ ] 无「复制粘贴重复逻辑」应复用而未复用的明显债务（同层内）

## E. 架构

- [ ] 改动落在正确 Layer 目录
- [ ] 无低层承载高层职责（如 L1 拼名称、L5 重选母体）
- [ ] 无连续同分下「只改 L5」糊名
- [ ] 无跨层错误 import / 语义越界

## F. 化学与硬约束

- [ ] 无 SMILES 特判 / 完整分子→名称表
- [ ] 无居所取代型投机
- [ ] 无 ML/LLM 命名；cache 未超 100
- [ ] 母体/位次/中英文组装无一眼违规（对照意图与 IUPAC 线索）
- [ ] 测试未放宽 normalize、未删断言换绿

## G. 裁决与共识

- [ ] 报告含：选题同意/意见、共识状态、代码质量、范围、Layer、门禁、裁决
- [ ] 裁决为 PASS | FIX | ROLLBACK | CHALLENGE 之一
- [ ] **PASS** 仅当选题同意 **且** 代码合格
- [ ] FIX 时「必须修改」可执行且未开新规则簇
- [ ] ROLLBACK 时写明 dual 或硬约束原因
- [ ] CHALLENGE 时含替代选题；未共识 **不** 建议 commit
- [ ] 往返 ≥3 仍无交集 → `escalation: human`
- [ ] PASS 时给出 workstate 可用摘要行

## 快速命令

```bash
python tools/structure_lint.py --root src/namepredict
pytest tests/unit -q
python -m benchmarks.benchmark --data data/merged_benchmark.json
```
