# chem-tdd-skill 检查清单

每轮 Agent **开始改代码前**与 **请求 commit 前**各过一遍。任一项未满足不得进入下一步。

## A. 开始前

- [ ] 已读本轮失败簇摘要与建议 Layer
- [ ] 已打开相关 `docs/iupac` 条款
- [ ] 已选择一个 **范围足够** 的改进点

## B. 抽例

- [ ] 正例 SMILES：**3–10** 个，来自失败簇
- [ ] 负例：**≥1** 个近邻（不应被本规则误伤）
- [ ] 期望 en/zh 来自benchmark原文或已验证正确输出；

## C. 测试（先红）

- [ ] 文件在 `tests/unit/`
- [ ] 文件头含 `# IUPAC: P-...` 与 `# Layer: Lx`
- [ ] 使用 `normalize_en` / `normalize_zh` 双语断言（中文可按 `zh is not None` 跳过）
- [ ] `assert r.success` 后再比名称
- [ ] 实现前运行：`pytest tests/unit -q` → **新增用例失败（红）**
- [ ] 未删除/弱化旧断言来制造假绿

## D. 实现（后绿）

- [ ] 仅修改允许路径：`src/namepredict/`（入口、layer0–5、cache）与 `tests/unit/`
- [ ] 无单 SMILES 特判（`if smiles == "..."` 等）
- [ ] 未改金标与 `benchmarks/benchmark.py`
- [ ] 单函数 ≤10 行、单文件 ≤1000 行；cache ≤100
- [ ] 实现后：`pytest tests/unit -q` → **全绿**

## E. 门禁命令

- [ ] `pytest tests/unit -q`
- [ ] `python tools/structure_lint.py`
- [ ] `python -m benchmarks.benchmark --data data/merged_benchmark.json`（或以编排器 `bench_cmd` 为准）
- [ ] dual 未下降 >0.5%；有提升（dual↑ 或 fails↓）才建议 commit

## F. 回合小结

- [ ] 已按 SKILL.md「回合小结格式」输出
- [ ] 含 IUPAC、Layer、TDD 红→绿、bench 前后、决策建议

## 快速命令

```bash
pytest tests/unit -q
python tools/structure_lint.py
python -m benchmarks.benchmark --data data/merged_benchmark.json
```
