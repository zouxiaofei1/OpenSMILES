---
description: Repo Wiki — 为 chem 项目自动生成结构化文档，支持增量更新、人工编辑保护、wiki_plan.yaml 前置干预。触发词：/repowiki、生成wiki、更新wiki、wiki文档、代码文档。
---

# Repo Wiki Skill

为 `namepredict`（SMILES → IUPAC 双语命名引擎）生成并维护结构化代码文档。

## 核心原则

1. **代码为源** — 所有文档声称必须有 `file:line` 源码锚点
2. **增量优先** — 更新时只改受影响页面，不全量重建
3. **人工优先** — 标注 `<!-- curated -->` 的段落永不被覆盖
4. **双语文档** — 中文为主，关键术语保留英文
5. **交叉引用** — 页面间用 `[[page-name]]` 链接
6. **最小约束** — 仅 `architecture/` 结构预定义；concepts、guides、reference 等目录由生成阶段自生长，不做预设

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
1. 读取 `wiki_plan.yaml`（如存在）
2. 列出该层所有 `.py` 文件，逐个阅读关键文件
3. 写作架构文档，必须包含：
   - **概述** — 该层在流水线中的位置、职责、输入/输出数据类型
   - **核心逻辑** — 关键函数/类的实现思路，≥500 字
   - **文件清单** — 该层所有源文件及职责说明
   - **数据流图** — 至少 1 个 Mermaid 图（flowchart/sequence/class）
   - **源码引用** — ≥3 处 `file:line` 锚点
   - **对外接口** — 该层 export 了哪些公开函数/类
   - **相关页面** — `[[wikilink]]` 到其他可能的页面（即使尚未创建）
4. 写完后返回该层发现的**跨层概念**列表（供 Phase 2 使用）

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

#### 质量门禁

每个 architecture 页面需满足：
- ≥500 字实质性内容（不含代码块/引用）
- ≥1 个 Mermaid 图
- ≥3 处源码锚点
- 明确说明该层与其他层的接口（谁调用它、它调用谁）

---

### `/repowiki:update` — 增量更新

1. 执行 `git diff --name-only HEAD` 获取变更文件
2. 与 `docs/wiki/log.md` 中记录的页面→源文件映射比对
3. 仅重写受影响的页面
4. 跳过所有 `<!-- curated -->` 标注的段落
5. 更新 `log.md`

---

### `/repowiki:modify <page>` — 修改指定页面

对话式修改，输入自然语言描述变更内容。保留所有 `<!-- curated -->` 段落。

---

### `/repowiki:supplement <page>` — 追加内容

向已有页面追加新内容，不删除现有内容。

---

### `/repowiki:rewrite <page>` — 完全重写

忽略所有保护标记，彻底重写指定页面。使用前需确认。

---

### `/repowiki:plan` — 创建/编辑 wiki_plan.yaml

生成 `wiki_plan.yaml` 模板或交互式编辑现有配置。

---

### `/repowiki:lint` — 健康检查

检查：孤立页面、失效链接（源码锚点指向不存在的文件或行）、过时页面（源文件更新时间 > 页面更新时间）。

---

## wiki_plan.yaml

位置：`docs/wiki/wiki_plan.yaml`

```yaml
version: 1

repowiki:
  # 预制模板: architecture | full_reference
  template: architecture

  # 生成阶段引导提示
  notes:
    - text: "重点关注 layer2 母体选择逻辑和 layer5 名称组装"
      author: "chem-dev"

  # 页面白名单（提供时严格按列表生成）
  # documents: []

  # 输出语言: zh | en
  language: zh

scope:
  include:
    - "src/namepredict/**"
  exclude:
    - "**/__pycache__/**"
    - "**/*.pyc"
```

关键配置：
- `notes` — 注入引导提示，让 AI 在生成时关注指定重点
- `documents` — 页面白名单，提供时**严格按列表生成**，不生成其他页面
- `scope` — 控制可见的文件范围（`.gitignore` 语法）

---


## 源码引用格式

```markdown
> **源:** `src/namepredict/layer2/parent_selector.py:142-158`
> 母体选择主逻辑，按官能团优先级遍历候选母体链。
```

---
### IUPAC 规则覆盖

项目实现了蓝皮书 P-1 至 P-49 的大部分规则，详见 `docs/RULES.md`。
