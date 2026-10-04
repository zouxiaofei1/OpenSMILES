---
description: Repo Wiki — 为 chem 项目自动生成结构化文档，支持增量更新、人工编辑保护、wiki_plan.yaml 前置干预。触发词：/repowiki、生成wiki、更新wiki、wiki文档、代码文档。
---

# Repo Wiki Skill

为 `opensmiles`（SMILES → IUPAC 双语命名引擎）生成并维护结构化代码文档。

## 文档结构（仅约束 architecture）

```
docs/wiki/
├── index.md                  # 总目录 + 项目概览（生成阶段最后写）
├── log.md                    # 操作日志（append-only）
├── wiki_plan.yaml            # 前置干预配置
└── architecture/             # ← 唯一预定义结构
    ├── overview.md           # 整体架构：6层流水线 + 数据流
    ├── layer0-preprocessor.md
    ├── layer1-analyzer.md
    ├── layer2-parent-selector.md
    ├── layer3-numbering.md
    ├── layer4-locants.md
    └── layer5-name-assembly.md
└── guides
└── reference
```

## 命令

以下命令均在对话中触发。

---

### `/repowiki:generate` — 从零生成 Wiki（subagent 并行）

**执行流程：两阶段 subagent 编排**

#### Phase 1: 并行写 architecture（6 subagent 同时跑）

为每一层启动一个 subagent，各自负责该层的架构文档：

```
layer0-agent  →  docs/wiki/architecture/layer0-preprocessor.md
layer1-agent  →  docs/wiki/architecture/layer1-analyzer.md
layer2-agent  →  docs/wiki/architecture/layer2-parent-selector.md
layer3-agent  →  docs/wiki/architecture/layer3-numbering.md
layer4-agent  →  docs/wiki/architecture/layer4-locants.md
layer5-agent  →  docs/wiki/architecture/layer5-name-assembly.md
```

每个 subagent 的任务：
1. 列出该层所有 `.py` 文件，逐个阅读关键文件
2. 写作架构文档
3. 写完后返回该层发现的**跨层概念**列表（供 Phase 2 使用）

#### Phase 2: 补充页面 + 写 index（协调 Agent）

Phase 1 全部完成后，启动 1 个协调 Agent：

1. 读取所有 6 份 architecture 页面
2. 汇总各 agent 报告的跨层概念
3. 自行决定需要哪些额外页面——例如：
   - 发现多处涉及官能团优先级 → 创建 `concepts/functional-group-priority.md`
   - 发现扩展点有固定模式 → 创建 `guides/adding-new-functional-group.md`
   - 发现核心数据类型被多层引用 → 创建 `reference/core-types.md`
   - **不要被目录名限制**——可以直接在 `docs/wiki/` 下创建任何有意义的页面
4. 为每个额外的页面启动 subagent 撰写（可并行）
5. 最后写 `index.md`（项目概览 + 完整页面地图）
6. 更新 `log.md`


### `/repowiki:update` — 增量更新

1. 执行 `git diff --name-only HEAD` 获取变更文件
2. 与 `docs/wiki/log.md` 中记录的页面→源文件映射比对
3. 仅重写受影响的页面
4. 跳过所有 `<!-- curated -->` 标注的段落
5. 更新 `log.md`

---


