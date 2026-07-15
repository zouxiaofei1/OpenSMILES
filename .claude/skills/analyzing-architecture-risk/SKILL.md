---
name: analyzing-architecture-risk
description: Use when asked whether a codebase will become hard to maintain, for architecture risk / tech-debt / maintainability audits, "难维护" "可维护性" "架构风险" "技术债" reviews, pre-refactor health checks, or before scaling features on a layered or agent-loop project. Triggers on symptoms like dual-track migrations, god layers, stringly IR, kind/enum explosion, docs drift, or half-wired abstractions.
---

# analyzing-architecture-risk

**Core principle:** Maintainability risk is measured from structure and process, not narrated from taste. No metric without a command or file cite; no verdict without a ranked matrix.

## When to use / not

**Use:** 用户问未来是否难维护、架构风险、技术债、扩展成本、是否该重构；功能工厂/分层流水线/双轨迁移项目的健康检查。

**Not:** 单文件风格点评；只要「这段函数怎么拆」；已明确只要实现计划且不要风险评估（可提示先做风险再规划）。

## Iron rules

1. **禁止空口指标** — 未跑命令/未读文件不得写 LOC、dual%、import 数、文件数。
2. **禁止只谈代码品味** — 必须覆盖：结构、契约、扩展轴、半迁移、文档、流程闸门。
3. **禁止只讲故事** — 交付物必须含风险矩阵 + 轨迹情景 + 按杠杆排序的建议。
4. **禁止把 lint 全绿当成低风险** — lint 管得了长度/import 方向，管不了组合爆炸与双轨。
5. **证据 > 资深意见** — 「别量 LOC」时仍要量；可在报告里注明意见与数据冲突。

## 强制流程（按序，可并行工具）

### Phase A — 地图（不写结论）

| 步 | 做什么 | 最低证据 |
|----|--------|----------|
| A1 | 入口与流水线 | 读主入口（如 `namer.py`/`main`），画出实际调用链 |
| A2 | 规模与分布 | 按模块/层统计 **文件数 + LOC**；标出顶格/最大文件 |
| A3 | 约束与门禁 | 读项目规则（prompt/CLAUDE/lint 脚本）：行数上限、函数上限、层规则、回归闸门 |
| A4 | 架构真源 | 对照 **代码树** vs **ARCHITECTURE/README/plans/workstate**；标漂移 |

### Phase B — 结构债（可工具化）

| 步 | 探测 | 风险信号 |
|----|------|----------|
| B1 | 跨层 import | 上层依赖下层 **私有** `_api`；下层 import 上层 |
| B2 | 扩展轴 | 新功能是否要求同时改 N 层 + 新枚举/kind/分支 |
| B3 | 契约形态 | dict string key / 无类型 IR vs dataclass/schema；key 数量级 |
| B4 | 半迁移 | 搜「零行为」「未接」「TODO wire」「legacy + new path」；新旧双轨 |
| B5 | 碎片化 | 私有微函数占比、同构复制模块（一实体一文件模式） |
| B6 | 测试与回归 | 测试是否按规则/簇组织；是否存在强 dual/coverage 闸门抑制重构 |

### Phase C — 过程债

| 步 | 探测 | 风险信号 |
|----|------|----------|
| C1 | 演进日志 | workstate/changelog：局部 +0.x% 是否主导、架构提交是否长期 dual 持平 |
| C2 | 闸门副作用 | 回退阈值/分数目标是否导致只敢「零行为架构」 |
| C3 | 知识载体 | 新人/新 session 会读到哪份文档？是否过期 |

### Phase D — 综合（必须写全）

报告 **必须** 含下列章节（可用中文标题）：

1. **现状快照表**（LOC/层分布/顶格文件/关键指标 + 来源）
2. **高 / 中 / 低风险**（每条：机制 → 症状 → 若不处理）
3. **风险矩阵**（严重度 × 紧迫度 × 症状 × 后果）
4. **轨迹情景** ≥2（继续现状 vs 收口架构；可选「过早铺全覆盖」）
5. **建议按杠杆排序**（P0/P1…，可执行，非空话）
6. **总判断**（现在难不难 / 未来会否更难 / 是否已注定 / 最危险时间点 / 最值得做的一件事）
7. **证据附录或「我做了什么」**（命令与结论对应；跳过的步骤显式写出）

## 通用风险模式速查（命中则升高评级）

| 模式 | 表现 | 典型后果 |
|------|------|----------|
| **组合爆炸扩展轴** | kind/type/enum × 骨架 × FG 各写 `_try_*` | 边际成本单调升 |
| **双轨半迁移** | 新 IR/engine 旁路，旧路径仍主跑 | 认知分裂、永久分叉 |
| **字符串 IR** | 跨层几十上百个 dict key | 重构无安全网 |
| **层目录≠层契约** | 跨层 import 私有函数 | 边界名存实亡 |
| **硬约束碎片化** | 极短函数 + 文件顶格 + 合法再拆文件 | lint 绿、概念散 |
| **胖中心层** | 单层占 LOC 大头且 import 中枢 | 单点失稳 |
| **回归闸门锁债** | 任何重构怕掉分 | 技术债锁定 |
| **文档漂移** | ARCHITECTURE 仍写旧文件名 | onboarding 误导 |
| **同构复制模块** | 每环/每实体一文件同一套 try | 抽象从未形成 |

## 评级启发式

- **极高：** 扩展轴仍是特例枚举 **且** 通用管线未接线；或双轨已超过一个迭代仍双写。
- **高：** 无类型跨层契约 / 私有 API 穿透 / 中心层 >50% LOC 仍在长。
- **中：** 文档漂移、测试结构尚可但闸门抑制重构、复制模式局部存在。
- **低：** 有自动化 lint、清晰入口、测试按行为组织、抽象在删旧而非只加新。

**总判断句式（必选）：** 短期是否可控 + 中长期在何种条件下滑向难维护 + 是否仍有收口窗口。

## 常见借口 → 回应

| 借口 | 现实 |
|------|------|
| 「风险意见不需要量 LOC」 | 无规模分布只能得到万能鸡汤 |
| 「lint 全绿所以结构健康」 | lint 不测组合爆炸与半迁移 |
| 「用户只要分析不要计划」 | 矩阵+轨迹+P0 建议是分析的一部分，不是另开项目 |
| 「dual/bench 太慢可跳过」 | 可跳过跑分，不可跳过 **读闸门规则与近期 dual 轨迹** |
| 「从经验已知分层系统都难维护」 | 必须用 **本仓库证据** 证明机制 |
| 「架构文档有了就够」 | 必须以代码树为准做漂移对照 |
| 「零行为架构提交说明在重构」 | 未接线 = 风险未降；要标「脚手架债」 |
| 「私有函数都是实现细节」 | 被外层 import 的 `_fn` 已是公开契约 |

## Red flags — 停下补做

- 报告没有 **一张** 含数据来源的快照表
- 只夸「函数很短/有分层」或只骂「代码丑」
- 未提扩展新功能要改几个点
- 未搜半迁移/双轨信号
- 未对比文档与代码
- 建议全是「加强文档/多写测试」而无结构杠杆
- 出现未测量的精确百分比

**任一命中：回到 Phase A–C 补证据，再写 Phase D。**

## 深度调节

| 用户意图 | 最低完成线 |
|----------|------------|
| 快速判断 | A1–A4 + B2 + B4 + 矩阵简化 + 总判断 |
| 详细分析（默认） | A–D 全做 |
| 要收口清单 | D 之后追加 30/90 天迁移顺序（仍先风险后计划） |

## 与其它 skill

- 收口实施计划 → `writing-plans` / `executing-plans`
- 具体 diff 审查 → 项目 code-review skill（如 `chem-code-review`）
- 本 skill **不** 替代实现；只产出风险结论与优先杠杆

## 最小命令模板（按需裁剪）

```bash
# 规模
find src -name "*.py" | wc -l
find src -name "*.py" -exec wc -l {} + | sort -n | tail -30

# 层/包分布：按目录汇总 LOC

# 跨层/私有依赖
rg -n "from package\\.layer|import package\\.layer" src
rg -n "from .* import _" src

# 半迁移
rg -n "零行为|未接|not wired|legacy|TODO|FIXME|deprecated" src docs workstate.md

# 契约松散度：.get("key") / ["key"] 频次

# 门禁
python tools/structure_lint.py --root src   # 若存在
```

将 `src`/包名换成实际树。Windows 可用同等 Glob/Grep/Read。

## 成功标准

- 他人能靠报告中的路径与命令复现主要数字
- 最高风险对应 **机制**（如何变难）而非 **观感**
- P0 建议若执行会直接削弱矩阵中的极高/高风险项
