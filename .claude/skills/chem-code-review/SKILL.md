---
name: chem-code-review
description: Use when reviewing namepredict code changes after an implementer subagent (prompt.txt step 4 / agent loop), before git commit, when dual may regress, or when challenging whether the chosen improvement topic is worth doing. Triggers on layer0–layer5 diffs, topic selection disputes, structure_lint, function length, layer purity, SMILES special-cases, residence-style shortcuts, Chinese locants, or benchmark rollback.
---

# chem-code-review

本 skill 约束 **检查/复核 Agent**（workflow 第二 subagent）：在实现之后、commit 之前，对 **选题** 与 **diff** 做项目特有审查，并决定 **通过 / 要求修复 / 回滚 / 反对选题**。

**编排位置**：`prompt.txt` 步骤 4–4.1（含与主 Agent 的选题共识环）；与 `skills/chem-tdd-skill` 配对——后者管「怎么写」，本 skill 管「该不该做 + 写完能不能留」。

## 核心原则

1. **可自动化的交给 lint**：`python tools/structure_lint.py` 已覆盖行数、函数体 ≤10、SMILES `==`、向上层 import、cache 条目上限。审查重点是 **lint 判不了的判断题**。
2. **代码合格 ≠ 选题合格**：即使 lint/pytest/架构全绿，仍须独立质疑「这一刀值不值得、粒度对不对」。例如只做「甲基」而通用单烷基骨架同等成本可覆盖 methyl/ethyl/propyl… → 应 **CHALLENGE**，不得因代码漂亮而 PASS。
3. **一次只放行一个改进点**：与 workflow「每次只处理一个问题」一致；发现无关大改 → 要求缩回。
4. **分数说话**：dual 回退 **>0.5%** → **必须回滚**该改动，不得「先 commit 再修」。
5. **系统命名 > 捷径**：居所取代型、单分子特判、只在 L5 糊字符串，一律拒。
6. **未达成选题共识不得 commit**：`CHALLENGE` 后进入主 Agent ↔ review 往返，直到双方对「做什么」一致（见「选题质疑与共识」）。

## 何时用 / 不用

**用**：实现 subagent 改完 `src/namepredict/**` 或 `tests/unit/**` 之后；commit 前；benchmark 可疑时；主 Agent 刚选定改进点、需要第二意见时（可只做选题审查、尚无 diff）。

**不用**：纯读 docs；只改 `agent_loop/` 编排；人类明确要求。

## 强制步骤（按序）

0. **选题审查（先于或并行于代码，且不可跳过）**：读 `workstate.md`「下一步候选」、本轮意图、失败簇/dual 预期；按「选题质疑清单」判断是否反对或要求改粒度。若明显不该做 → 可直接 `CHALLENGE`，**不必**为凑报告强跑完整 bench（硬红旗/已实现 diff 仍建议跑门禁以决定是否 ROLLBACK 已改代码）。
1. **读 diff 范围**：`git diff` / 改动文件列表；确认仅触及允许路径。
2. **跑 structure_lint**（若实现方未跑或结果未给出）：
   ```bash
   python tools/structure_lint.py --root src/namepredict
   ```
3. **静态审查**（下文清单 + 层职责表）；逐函数量体、对层边界与化学逻辑做判断。
4. **跑 pytest**（若有 `tests/unit` 改动或实现声称 TDD）：
   ```bash
   pytest tests/unit -q
   ```
5. **跑 benchmark**（门禁轮次不可省；纯选题 `CHALLENGE` 且工作区无实现 diff 时可标「未跑」）：
   ```bash
   python -m benchmarks.benchmark --data data/merged_benchmark.json
   ```
   记录 **中英文双重准确: [ ]%**（dual）；与改前对比。
6. **裁决**：输出「审查报告」；`FIX`/`ROLLBACK`/`CHALLENGE` 时进入对应闭环，**禁止**代为实现无关新规则。

## 选题质疑与共识

审查 Agent **默认带怀疑**：选题由主 Agent 提出，不享有免审特权。代码质量与选题质量 **解耦评分**。

### 选题质疑清单（命中 ≥1 条应倾向 CHALLENGE）。选题应该足够有影响性，过小应当 CHALLENGE

| 信号 | 例 | 更优方向 |
|------|----|----------|
| **过窄切片** | 只做 methyl，实现已是「数碳+词干」却不推广到 C_n 烷基 | 通用 monoalkyl（methyl…）同一规则一次做完 |
| **过碎 / 违反自然规则边界** | 把同一 IUPAC 条款拆成多轮刷 commit | 一轮吃透一条主规则 |
| **预期 dual 增益过低** | 无结构推进 | 换失败簇更大或卡主路径的规则 |
| **层优先级错误** | 同分根因在 L2/L3，却选题 L5 词表 | 先结构性层 |
| **依赖未就绪** | 取代基命名但 L2 母体仍错，测例会假绿或无意义 | 先堵母体/编号 |
| **与 workstate 候选无对话** | 忽略已写「下一步候选」且无反驳 | 采纳候选或书面说明为何改道 |
| **特判气质选题** | 「先让这 3 个 SMILES 过」 | 拒绝；改系统规则选题 |
| **范围膨胀** | 名义一个点，diff 横跨多官能团多条款 | 缩回单主规则 |



### 反对时必须给出的内容

`CHALLENGE` 时给出：

1. **反对什么**：当前选题一句话  
2. **为何不值 / 为何粒度错**：对应上表信号 + 与 dual/架构的关系  
3. **替代选题（1–2 个）**：更合适的主规则与 Layer，仍单点  
4. **若坚持原选题的条件**：例如「扩成 C1–C4 直链烷基同一实现才接受」  
5. **已有代码怎么处理**：保留改 / 要求改范围重写 / 工作区回滚后再开题  

### 共识环（直到一致）

```text
主 Agent 选题（或已实现）
    → review：选题 + 代码审查 → 报告
         ├─ PASS（选题同意 + 代码合格）→ 主 Agent 可 commit
         ├─ FIX → 主 Agent 按「必须修改」改代码 → 再 review（可只复检 FIX 项）
         ├─ ROLLBACK → 回滚工作区 → 主 Agent 重新选题或换方案 → 再 review
         └─ CHALLENGE → 主 Agent 书面回应（接受替代 / 辩护原题 / 折中粒度）
                              → review 再裁决（同意则撤 CHALLENGE 改 PASS/FIX；
                                 仍不同意则继续 CHALLENGE，理由须新增或收紧条件）
                              → 循环直到选题一致；**一致前禁止 git 提交**
```

**主 Agent 辩护原题时**须答：预期 dual/fails 影响、为何不能通用化、与 workstate 关系、S3 边界。空泛「先做简单的」→ review 可继续反对。

**僵持**：往返 ≥3 轮仍无交集 → 报告标 `CHALLENGE` + `escalation: human`，暂停交给人类，**仍不得 commit**。

**角色边界**：review 可反对选题、要求换题或扩/缩范围；**不**自己开写新规则实现。主 Agent 负责改代码或换题。

## 可写 / 只读（复核时再确认）

| 允许出现在 diff | 禁止出现在 diff |
|----------------|----------------|
| `src/namepredict/` 入口、`layer0`–`layer5`、受限 `cache` | `data/*`、`data/merged_benchmark.json` |
| `tests/unit/` | `benchmarks/benchmark.py` 及计分逻辑 |
| （编排需要时）`workstate.md` / `agent_loop/memory/` | 其它 layer 文件混入错误 layer 目录 |

`namepredict` 下 **非** 入口与 layer0–5 / 约定 cache 的文件视为只读。临时文件只应在 `tools/`（或项目约定的 chem/tools），**不应**进生产包。

## Layer 职责（越界 = 架构问题）

| Layer | 应做 | 不应做 |
|-------|------|--------|
| **L0** Preprocessor | RDKit 清洗、H、盐、无机/有机分流 → 分子单元 | 命名、选母体、写 en/zh 字符串 |
| **L1** Analyzer | 环/链/官能团检测 → `MolecularGraph` 类信息（rings, chains, fg_list…） | 选最优母体、位次、拼名称 |
| **L2** Parent Selector | 按 IUPAC P-44/P-45 评分选母体 → `ParentContext` | 拼最终名；把取代基完整命名逻辑塞进来 |
| **L3** Substituent Extractor | 非母体片段类型、连接点、递归取代基、优先级 | 全局编号方案；最终 en/zh 组装 |
| **L4** Numbering | 候选编号、最低位次集、字母序次级、立体影响 → `NumberingResult` | 中英文词干表；cache 俗名 |
| **L5** Assembler | 按编号结果组装 en/zh（位次、括号、后缀） | 重新选母体/重算最长链；SMILES 特判 |

**硬边界**：

- 不得在 `layerN/` 目录放入属于其它 layer 的职责文件。
- 不得 `layerK` import `layerM`（M>K）——lint 已抓显式 import；审查时还要抓 **把高层逻辑复制粘贴进低层**。
- 连续同分失败却只改 L5 字符串 → **架构性红旗**，要求回到 L1/L2/L4。

## 函数式与体量（lint 之外）

`structure_lint` 数的是函数体 **非空非注释行 ≤10**。审查额外要求：

| 检查项 | 通过标准 |
|--------|----------|
| 单函数职责 | 一个明显意图；复杂控制流应拆成 `_helper` |
| 输入输出 | 尽量单一主输入/主输出；避免巨型 kwargs 字典无文档 |
| 复用 | 与邻层/同文件已有 `_foo` 重复的逻辑应复用或下沉共享小函数（仍须落在正确 layer） |
| 单文件 | ≤500 行；逼近时要求拆文件但 **不跨层** |
| 改动后超限 | **必须拆分**后再通过，不得「先过门禁以后再拆」 |

**拆分配方**（审查时要求实现方按此改，而不是你扩写规则）：

```text
过长函数 → 按「纯计算步」抽出 _step_a / _step_b
          → 入口函数只编排调用（仍 ≤10 行）
          → 每个 _helper 同样 ≤10 行
```

**禁止实现类**（发现即 **阻断**）：深度学习/LLM 命名；改 benchmark 数据库或计分；cache 无界增长。

## 回归门禁

| 门禁 | 命令/来源 | 失败动作 |
|------|-----------|----------|
| structure_lint | `python tools/structure_lint.py --root src/namepredict` | 不通过 → 修复前不得 commit |
| pytest | `pytest tests/unit -q` | 红 → 阻断 |
| benchmark dual | `python -m benchmarks.benchmark --data data/merged_benchmark.json` | dual **下降 >0.5%** → **回滚** |
| 提升判定 | dual↑ 或 fails↓ | 无提升且无修复价值 → 不建议 commit |

报告分数时使用：**中英文双重准确: X%**（与 `prompt.txt` / workstate 一致）。对比必须用 **同一 bench 命令** 的改前/改后。

## 输出：审查报告（固定形状）

审查结束时 **只输出** 下列结构（可中英混排术语，正文简体中文）：

```markdown
## 代码审查报告
- **范围**: <改动文件相对路径列表；无 diff 则写「仅选题」>
- **主规则/意图**: <一句话；对应 IUPAC 条款若可知>
- **选题同意**: 是 / 否
- **选题意见**: <同意理由；或反对信号 + 替代选题 1–2 个 + 坚持原题条件>
- **共识状态**: 一致 / 待主 Agent 回应 / 僵持→human（往返轮次: n）
- **Layer**: <L0–L5；是否越界>
- **structure_lint**: pass/fail/未跑（摘要）
- **pytest**: pass/fail / 未跑+原因
- **benchmark**: dual <before>% → <after>%；**中英文双重准确: <after>%** / 未跑+原因
- **体量**: 函数体超限 <n> 处；文件超限 <n> 处（明细路径:行）
- **架构**: 无问题 / 问题列表（层混用、高层逻辑下沉、L5 糊名…）
- **化学逻辑**: 无问题 / 问题列表（居所捷径、特判、位次、母体链…）
- **路径合规**: 可写范围内 / 违规路径 / 无 diff
- **代码质量**: 合格 / 不合格（与选题解耦；可「代码合格但选题否」）
- **裁决**: PASS | FIX | ROLLBACK | CHALLENGE
- **必须修改**（FIX 时逐条可执行）:
  1. ...
- **回滚原因**（ROLLBACK 时）: dual 回退 ... / 硬约束违反 ...
- **挑战要点**（CHALLENGE 时）: 反对点；替代选题；对已有 diff 的处置建议
- **可提交摘要**（仅 PASS 时，供 workstate）:
  [#pending][IUPAC ...] 标题 [+n tests, dual a%→b%]
```

**裁决定义**（可组合理解，最终只选一个主裁决）：

| 裁决 | 何时 |
|------|------|
| **PASS** | **选题同意** 且 代码合格：lint/pytest 绿；dual 未回退 >0.5%；无架构/硬约束问题；体量与路径合规 |
| **FIX** | 选题同意（或选题争议已关闭），但代码可局部改正；给出有序列表，**不要**顺手开新规则 |
| **ROLLBACK** | dual 回退 >0.5%，或改金标/特判/ML/越权路径等硬违反（与选题是否聪明无关） |
| **CHALLENGE** | 反对选题或粒度（**即使代码完全合格**）；或主 Agent 辩护不成立；一致前禁止 commit |

优先级：存在硬红旗/ dual 崩 → **ROLLBACK** 优于 CHALLENGE；选题不同意且代码仅因错题而存在 → **CHALLENGE**（可附带「建议回滚工作区」）；选题同意仅代码差 → **FIX**。

## 常见借口（审查时直接驳回）

| 借口 | 现实 |
|------|------|
| 「lint 过了就行」 | lint 不看层语义、居所捷径、中文位次、dual、选题 |
| 「代码没问题为什么不 PASS」 | 代码合格 ≠ 选题合格；过窄/错层选题 → CHALLENGE |
| 「只做甲基，下一步再做乙基」 | 若实现已是通用烷基骨架，拆轮是刷 commit；要求一次规则化 |
| 「简单点先落地」 | 须证明是最小可证伪步，而非回避自然规则边界 |
| 「只有 12 行，差不多 10」 | ≤10 是硬顶；超了就拆 |
| 「先 commit 分数掉了再修」 | >0.5% 必须当轮回滚 |
| 「这是临时特判，以后系统化」 | 特判进主路径即债务；cache 仅俗名且有上限 |
| 「只改了 L5，最快」 | 同分结构性问题禁止只补 L5 |
| 「测试太严，放宽 normalize」 | 禁止；双语严格断言 |
| 「顺手把旁边失败也修了」 | 违反单改进点 / S3；要求缩回 |
| 「benchmark 太慢，跳过」 | 有 diff 的门禁轮次须跑 bench；跳过不得 PASS |
| 「review 不能否选题」 | 可以；共识前不得 commit |

## 红旗 — 立即 STOP

- diff 触及 `data/*` 或 `benchmarks/benchmark.py`
- `if smiles ==` / 完整 SMILES→名称表（生产路径）
- 函数体 >10 仍声称完成
- dual 下降 >0.5% 仍建议 commit
- layer 目录职责明显错位或跨层文件
- 居所取代型或「枚举分子出名称」
- 为通过而删测试、改金标、放宽匹配
- **选题未一致仍准备 git commit**
- **空泛辩护后再次原题提交且无新增理由**

**任一红旗 → 不得 PASS。**

## 与 chem-tdd-skill 的分工

| | chem-tdd-skill | chem-code-review（本 skill） |
|--|----------------|------------------------------|
| 角色 | 实现 Agent | 检查 Agent |
| 重点 | 先红后绿、抽例、S3、回合小结 | **选题质疑**、层边界、体量拆分、捷径、dual 回滚 |
| 选题 | 主 Agent 提出并实现 | 可反对；替代选题；直到共识 |
| 测试 | 必须自己写并证明红→绿 | 核实仍在且未放宽；不替实现方新开规则簇 |
| bench | 实现后跑 | **再跑**并执行 >0.5% 回滚（纯选题挑战且无 diff 可免） |

## 快速检查清单

- [ ] 选题已审；过窄/错层/低收益等已表态（同意或 CHALLENGE）
- [ ] 路径仅允许区；无金标/计分改动
- [ ] structure_lint 绿（有 diff 时）
- [ ] 每函数体 ≤10；超限已要求拆分
- [ ] 层职责表无越界；无「只 L5 糊同分」
- [ ] 无 SMILES 特判 / 无居所捷径 / cache 未爆
- [ ] 中文位次与 en 组装不明显违规
- [ ] pytest 绿（若适用）
- [ ] dual 对比完成；回退 ≤0.5%（有 diff 时）
- [ ] 裁决与报告格式完整；PASS 仅当选题同意且代码合格
- [ ] CHALLENGE 含替代选题与 diff 处置；未共识不 commit
