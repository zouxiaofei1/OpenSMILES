# 新增环系指南

本文档分步说明如何为 NamePredict 流水线新增一个环系（ring system）——无论是新的 fused heterocycle（稠杂环）、bridged system（桥环）、spiro system（螺环），还是新的单杂环。

---

## 前置知识

阅读本指南前，建议先了解以下架构文档：

- [[architecture/overview]] — 六层流水线总览
- [[architecture/layer2-parent-selector]] — Layer2 母体选择器的核心架构（候选生成 + 评分排序）
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 候选管线，无需注册 orienter）
- [[architecture/layer5-name-assembly]] — Layer5 中英双语名称组装
- [[concepts/functional-group-priority]] — FG 优先级规则

---

## 流程概览

新增一个环系涉及以下改动，按顺序为：

1. **Layer1 环检测**（通常无需改动，SSSR 已自动识别）
2. **Layer2 骨架识别**（主要工作量——在 `ring_scaffold.py` 注册 ScaffoldSpec / 保留拓扑）
3. **Layer4 编号**（**通常无需改动**——`numbering_engine` 的 P-14.4 候选管线自动适用）
4. **Layer5 名称组装**（添加命名函数 / `_KIND_TABLE` 词干）
5. **FG-环组合支持**（如环上有 COOH/CHO/CN/OH/NH2 等 FG 变体）

> 核心原则（2026-08 正交化）：**scaffold × kind × 数量全正交**——加新环只加 `ScaffoldSpec` + 词干表，不再枚举 cycloalcohol/naphthalenol/indolol 组合 kind。

---

## 第一步: 环检测前提（Layer1）

大多数环系已被 Layer1 的 `ring_systems.py` 通过 **SSSR（Smallest Set of Smallest Rings）+ Union-Find 融合** 自动检测，无需额外工作。

> **源:** `src/namepredict/layer1/ring_systems.py`

SSSR 利用 RDKit 的 `GetRingInfo().AtomRings()` 获取所有最小环，然后通过 fusion graph（共享 >=2 原子的边）将稠合环合并为 ring system。桥环（共享 >=3 原子）和螺环（共享 1 原子）也被识别。

如果你的环系需要**超越 SSSR 的特殊识别**，则需扩展 ring_systems.py 或创建新的检测模块。

---

## 第二步: 母体选择（Layer2）—— 主要工作量

### 2.1 骨架识别入口

母体选择统一走 P-44 规则管线（`rule_driven_parent_candidates`），环骨架身份由 **`ring_scaffold.py` 的 `resolve_ring_scaffold`**（`:287-294`）在表达阶段解析：

```
resolve_ring_scaffold(info, skeleton)
  ├─ 1. get_identity(skeleton.scaffold_id)  # _ALL_SPECS 直接命中
  ├─ 2. _matched_id → _TOPOLOGY.match_systems  # 按 ring_system 原子集合匹配保留拓扑
  └─ 3. _generic_carbocycle  # 全碳非保留环 → ScaffoldIdentity("carbocycle",...)
```

### 2.2 注册 ScaffoldSpec（新的命名骨架）

如果环系是 retained-name ring（如吡啶、萘、吲哚）或新的 scaffold，需在 **`src/namepredict/layer2/ring_scaffold.py`** 的 `_ALL_SPECS` 中添加一个 `ScaffoldSpec`。

`ScaffoldSpec` 字段（`ring_scaffold.py:34-50`）：

| 字段 | 类型 | 说明 |
|------|------|------|
| `id` | `str` | 与 parent dict 的 `kind` 一致 |
| `naming_class` | `str` | 命名类别，如 `"monohetero"`, `"fused56"`, `"carbocycle"`, `"naph_family"` |
| `stem_en` / `stem_zh` | `str` | 环系骨架的英文/中文 stem |
| `n_rings` | `int` | 环数（单环=1, 双环稠合=2, 三环=3） |
| `ring` | `str` | `"hetero"` 或 `"carbo"` |
| `retained` | `bool` | 是否为 IUPAC 保留名 |
| `fg_rank` | `int` | 官能团优先级（0 表示纯母体） |
| `numbering` | `NumberingPolicy` | 编号策略 |
| `sub_rules` / `principal_slots` | — | 取代基/主官能团插槽（可选） |

`ScaffoldSpec` 会被 `kind_registry.py` 的 `_load_from_scaffold_specs`（`:141`）自动加载为 `KindMeta`（bootstrap 最后一步，覆盖同 kind 已有注册），提供 scoring 和 stem 元数据。

> **源:** `src/namepredict/layer2/ring_scaffold.py`, `src/namepredict/layer2/kind_registry.py:141`

### 2.3 注册保留母体（唯一事实来源 `_TEMPLATES`）

环骨架身份由 `resolve_ring_scaffold` 识别：① `get_identity(skeleton.scaffold_id)` 显式命中 → ② **SMILES 模板子图同构**（`match_retained`）→ ③ 兜底 `_generic_carbocycle`（全碳非保留环 → carbocycle）。

新增保留环系（如 quinoline）**只改 `ring_scaffold.py` 的 `_TEMPLATES` 一张表**——加一条 `{smiles, stem_en, stem_zh, naming_class}`。ScaffoldSpec 由 `_spec_from_template` 自动派生（n_rings/ring 从 smiles 算，retained=True，fg_rank=0），`kind_registry` 经 `all_specs()` 自动注册词干。位置异构体（quinoline/isoquinoline、二嗪、二唑等）在元素标注的子图同构下天然区分，无需额外消解。

> `_TOPOLOGY` 五元组表与手写 `_ALL_SPECS` 已删除；`match_systems`/`match_scaffold_ids`/`registry`/`get_entry` 为模板语义查询。

### 2.4 FG-环组合的 kind（无需手动注册组合 kind）

**不再需要**为 FG-环组合（如 pyridinol、cycloalkanol、naphthalenol）手动注册组合 kind——2026-08 正交化已根除组合 kind。环上带 FG 时：

- **环 + 主 FG**：`express_ring_principal` 把 kind **收敛为 FG 类别**（alcohol/acid/amine/...），词干由 scaffold 承载（`pack_parent_stem` 按 scaffold_id 注入）；苯/饱和环/稠环/杂环一律平等
- **苯 + FG 保留名**：chain_engine variant 提供 phenol/benzoic acid/aniline/benzaldehyde/benzonitrile/benzamide 等（苯专属）；其余 scaffold 走通用词干命名（naphthalen-1-ol / pyridine-3-carboxylic acid / imidazol-2-amine）
- 命名时 `chain_engine._KIND_TABLE` 的 entry 按 `scaffold_id` 注入词干（`_ring_stem` 从 parent 的 stem_en/zh 派生，IUPAC 词干去尾部 e，苯除外）

---

## 第三步: 编号（Layer4）—— 通常无需改动

**`numbering_engine.orient_numbering`**（P-14.4 候选管线）不再按 kind 枚举 orienter：

1. 保留 scaffold 经 `plan_from_chain` 用 L2 注入的 `numbering_scaffold` 事实直接定向
2. 其余环枚举 2n 个候选（每原子 1 号位 × 双向），按 principal FG → 多重键 → 取代基位次集逐条收窄

因此**新增环系无需写 orienter**。唯一需要保证的是 L2 正确注入 `scaffold_id` 与（对保留 scaffold）`numbering_scaffold` 事实。如果环系需要固定编号（如稠环的 `3a`/`4a` 标签），在 `numbering_scaffold_facts`（`ring_scaffold.py:147-155`）提供对齐的 `{scaffold_id, labels}` 事实即可。

> **源:** `src/namepredict/layer4/numbering_engine.py`, `src/namepredict/layer4/locants/adapt.py`

---

## 第四步: 名称组装（Layer5）

### 4.1 保留名环系

对于保留名环系（苯、吡啶、萘、吲哚），词干在 `ring_scaffold.py` 的 `ScaffoldSpec.stem_en`/`stem_zh` 定义。`chain_engine._KIND_TABLE` 的 entry 在 `_names_for` 中按 `scaffold_id` 运行时替换：

```python
# assembler.py: _RING_STEM = {naphthalene: ("naphthalen","萘"), indole: ..., pyridine: ..., quinoline: ...}
sid in _RING_STEM → replace(entry, stem=_RING_STEM[sid], ...)
sid == "carbocycle" → replace(entry, cyclic=True, ene_loc_omit=True, ...)
```

### 4.2 特殊命名需求

如果你的环系需要特殊的命名处理，则可能需要：

- 在 `layer5/benzene_names.py` 添加命名函数
- 在 `assembler_prefixes.py` 中注册前缀构建规则
- 在 `layer5/stems.py` 中添加 stem 表

---

## 第五步: FG-环组合支持

正交化后组合 kind 已根除，FG-环组合的命名由以下路径承担：

- **苯系保留名**（benzoic/phenol/aniline/benzaldehyde/benzonitrile/benzamide/benzoate）：`layer5/typed_kinds.py` 的 `_BENZENE_RETAINED`（要求 `scaffold_id=="benzene"` 且 `multiplicity==1`）
- **稠环/杂环 FG 收敛**：`_RING_FG_SCAFFOLDS = {naphthalene, indole, pyridine, quinoline}` → 收敛为 FG 类别，词干由 chain_engine 注入
- **环外酸**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_acid_names` → `cyclohexanecarboxylic acid`
- **环二酸立体化学**：`layer4/cyclo_relative_stereo.py`（`acid` + `scaffold_id=="carbocycle"` → cis/trans）

---

## 文件改动清单

| Layer | 文件 | 改动内容 |
|-------|------|---------|
| L1 | `layer1/ring_systems.py` | 特殊环检测（通常无需改动） |
| L2 | `layer2/ring_scaffold.py` | 新增 `ScaffoldSpec`（`_ALL_SPECS`）+ 保留拓扑条目（`_TOPOLOGY`） |
| L2 | `layer2/kind_registry.py` | 通常无需改动（`_load_from_scaffold_specs` 自动注册） |
| L4 | `layer4/ring_scaffold.py` (`numbering_scaffold_facts`) | 固定编号事实（如稠环 `3a`/`4a` 标签） |
| L5 | `layer5/assembler.py` | `_RING_STEM` 词干表 / worker |
| L5 | `layer5/benzene_names.py` 或 `layer5/chain_engine.py` | 命名函数 / `_KIND_TABLE` spec |
| L5 | `layer5/typed_kinds.py` | 环系 FG 收敛（`_RING_FG_SCAFFOLDS`） |
| L5 | `layer5/assembler_prefixes.py` | 前缀构建规则（如需） |
| — | `tests/` | 添加 SMILES 测试用例 |

---

## 工作示例: 新增吡啶（pyridine）支持

以下是 NamePredict 中吡啶支持的实现路径（实际代码），理解它即可掌握新增保留名单杂环的完整流程。

### Step 1 — 环检测

吡啶的六元芳环（5C + 1N）被 SSSR 自动检测。未对 Layer1 做额外改动。

### Step 2 — Layer2 骨架识别

pyridine 走**数据驱动**识别，全部在 `ring_scaffold.py`：

1. **ScaffoldSpec**：`_ALL_SPECS` 已有 `pyridine`（`naming_class="monohetero"`, `stem_en="pyridine"`, `stem_zh="吡啶"`, `retained=True`）
2. **保留拓扑匹配**：`_TOPOLOGY` 已有 `("pyridine", 1, 6, (7,), "mono", True)`，`resolve_ring_scaffold` 经 `match_systems` 标注 `scaffold_id="pyridine"`

### Step 3 — 编号

无需 orienter。`numbering_engine` 的 P-14.4 候选管线对环枚举 2n 个候选；`_ring_hetero_start`（`numbering_engine.py:119-125`）把单杂环的唯一杂原子固定为 1 号位（P-14.4，pyridine 的 N=1）。

### Step 4 — 命名

`assembler._names_for` 查 `_KIND_TABLE`（`alkane` entry），按 `scaffold_id=="pyridine"` 从 `_RING_STEM` 注入词干 `("pyridine","吡啶")`。FG 变体（pyridinecarboxylic 等）经 `typed_kinds` 收敛为 FG 类别后同样走 chain_engine。

### Step 5 — 测试

SMILES 测试用例示例：
- `c1ccncc1` — pyridine / 吡啶
- `c1ccncc1C(=O)O` — pyridine-2-carboxylic acid（吡啶甲酸 / 吡啶-2-甲酸）
- `c1ccncc1N` — pyridin-2-amine（吡啶-2-胺）

---

## 常见问题

### Q: 有 ScaffoldSpec 但没有骨架识别路径？

`ScaffoldSpec` 仅提供 stem 和编号策略的元数据；分子仍须在 `resolve_ring_scaffold` 的识别路径（`_ALL_SPECS` 直查 / `_TOPOLOGY` 匹配）中被命中，其骨架才会使用该 scaffold。

### Q: 为什么有些环系没有独立的 .py 文件？

部分简单环系（如 furan, thiophene, pyrrole）不设独立模块——它们经 `ring_scaffold.py` 的 `_ALL_SPECS` 批量注册，共享相同的编号策略，仅杂原子与 stem 不同。

### Q: 多环稠合系统如何处理？

多环稠合系统（如 naphthalene, indole）经 `ring_systems.py` 的 ring_system 分组识别，骨架身份由 `resolve_ring_scaffold` 解析（`_TOPOLOGY` 匹配）。编号方面，L2 注入 `numbering_scaffold` 事实（`{scaffold_id, labels}`），Layer4 的 `plan_from_chain` 直接定向。

### Q: ScaffoldSpec 和 KindMeta 的关系？

`ScaffoldSpec` 是 stem 和编号的单一权威来源（Single Authority）。`kind_registry.py` 在 bootstrap 的最后一步（`_load_from_scaffold_specs`，`:141`）读取所有 `ScaffoldSpec` 并转换为 `KindMeta` 注册。如果一个 kind 同时有 ScaffoldSpec 和手动 KindMeta 注册，ScaffoldSpec 覆盖前者（因为最后执行）。

> **源:** `src/namepredict/layer2/kind_registry.py:141`

---

## 相关文档

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器完整架构
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 候选管线）
- [[architecture/layer5-name-assembly]] — Layer5 名称组装
- [[concepts/functional-group-priority]] — FG 优先级与 IUPAC P-44 规则
