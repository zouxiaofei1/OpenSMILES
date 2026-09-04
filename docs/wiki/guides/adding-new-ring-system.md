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
2. **Layer2 骨架识别**（主要工作量——在 `ring_scaffold.py` 的 `_TEMPLATES` 注册模板）
3. **Layer4 编号**（**通常无需改动**——`numbering_engine` 的 P-14.4 候选管线自动适用）
4. **Layer5 名称组装**（添加 `_KIND_TABLE` 词干 / variant）
5. **FG-环组合支持**（如环上有 COOH/CHO/CN/OH/NH2 等 FG 变体）

> 核心原则：**scaffold × kind × 数量全正交**——加新环只加 `_TEMPLATES` 模板 + 词干表，不枚举 cycloalcohol/naphthalenol/indolol 组合 kind。

---

## 第一步: 环检测前提（Layer1）

大多数环系已被 Layer1 的 `ring_systems.py` 通过 **SSSR（Smallest Set of Smallest Rings）+ Union-Find 融合** 自动检测，无需额外工作。

> **源:** `src/namepredict/layer1/ring_systems.py`

SSSR 利用 RDKit 的 `GetRingInfo().AtomRings()` 获取所有最小环，然后通过 fusion graph（共享 >=2 原子的边）将稠合环合并为 ring system。桥环（共享 >=3 原子）和螺环（共享 1 原子）也被识别。

如果你的环系需要**超越 SSSR 的特殊识别**，则需扩展 ring_systems.py 或创建新的检测模块。

---

## 第二步: 母体选择（Layer2）—— 主要工作量

### 2.1 骨架识别入口

母体选择统一走 P-44 规则管线（`rule_driven_parent_candidates`），环骨架身份由 **`ring_scaffold.py` 的 `resolve_ring_scaffold`**（`:398`）在表达阶段解析：

```
resolve_ring_scaffold(info, skeleton)
  ├─ 1. get_identity(skeleton.scaffold_id)  # 模板 id 直接命中
  ├─ 2. match_retained(info, atom_ids)      # SMILES 模板子图同构精确覆盖
  └─ 3. _generic_carbocycle                 # 全碳非保留环 → ScaffoldIdentity("carbocycle",...)
```

### 2.2 注册保留母体（唯一事实来源 `_TEMPLATES`）

环骨架身份由 `resolve_ring_scaffold` 识别：① `get_identity(skeleton.scaffold_id)` 显式命中 → ② **SMILES 模板子图同构**（`match_retained`）→ ③ 兜底 `_generic_carbocycle`（全碳非保留环 → carbocycle）。

**新增保留环系（如 quinoline）只改 `ring_scaffold.py` 的 `_TEMPLATES` 一张表**——加一条 `{smiles, stem_en, stem_zh, naming_class}`：

```python
# ring_scaffold.py `_TEMPLATES`（`81` 起）
"quinoline": {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline", "stem_zh": "喹啉", "naming_class": "fused56"},
```

`_spec_from_template`（`:212`）自动派生 `ScaffoldSpec`（n_rings/ring 从 smiles 算，retained=True），`all_specs()`/`get_spec()`/`get_identity()` 均由此派生；`kind_registry._load_from_scaffold_specs`（`kind_registry.py:101`）据此自动注册 KindMeta 词干（bootstrap 唯一一步，ScaffoldSpec 是词干权威）。位置异构体（quinoline/isoquinoline、二嗪、二唑等）在元素标注的子图同构下天然区分，无需额外消解。

> 无 `_TOPOLOGY` 五元组表与手写 `_ALL_SPECS`；`match_systems`/`match_scaffold_ids`/`registry`/`get_entry` 为模板语义查询。

### 2.3 FG-环组合的 kind（无需手动注册组合 kind）

**无需**为 FG-环组合（如 pyridinol、cycloalkanol、naphthalenol）手动注册组合 kind——scaffold × FG 已正交化。环上带 FG 时：

- **环 + 主 FG**：`express_ring_principal` 把 kind **收敛为 FG 类别**（alcohol/acid/amine/...），词干由 scaffold 承载（`pack_parent_stem` 按 scaffold_id 注入）；苯/饱和环/稠环/杂环一律平等
- **苯 + FG 保留名**：chain_engine `_KIND_TABLE` 各 entry 的 `variant` 提供 phenol/benzoic acid/aniline/benzaldehyde/benzonitrile/benzamide 等（苯专属）；其余 scaffold 走通用词干命名（naphthalen-1-ol / pyridine-3-carboxylic acid / imidazol-2-amine）
- 命名时 `chain_engine._KIND_TABLE` 的 entry 按 `scaffold_id` 注入词干（`assembler._ring_stem` 从 parent 的 stem_en/zh 取完整词干；结尾 e 的省略由 chain_engine `_elide_parent_e` 按后缀首字母判定，P-60.2(a)——diol/dione/carbaldehyde 等辅音开头后缀保留 e）

---

## 第三步: 编号（Layer4）—— 通常无需改动

**`numbering_engine.orient_numbering`**（P-14.4 候选管线）不再按 kind 枚举 orienter：

1. 枚举候选编号：链正反（2 个）/ 环每原子 1 号位 × 双向（2n 个）
2. `_fixed_start` 固定 1 号位：杂原子环优先杂原子（Z 最小 = locant 1），否则 FG 锚点/自由基字段
3. 按 principal FG → 多重键 → 取代基位次集逐条收窄

因此**新增环系无需写 orienter**。唯一需要保证的是 L2 正确注入 `scaffold_id`。固定编号事实（稠环 `3a`/`4a` 标签）由 `numbering_scaffold_facts`（`ring_scaffold.py:254`）基于 `_TEMPLATES` 生成，L4 校验其存在性（`numbering_scaffold_required`）。

> **源:** `src/namepredict/layer4/numbering_engine.py`, `src/namepredict/layer2/ring_scaffold.py:254`

> 注：无 `NumberingPlan` 机制与 `locants/` 子包——编号完全走候选枚举 + `chain.index + 1`，不依赖固定编号 plan。

---

## 第四步: 名称组装（Layer5）

### 4.1 保留名环系

对于保留名环系（苯、吡啶、萘、吲哚等），词干在 `ring_scaffold.py` 的 `_TEMPLATES`（`stem_en`/`stem_zh`）定义。`chain_engine._KIND_TABLE` 的 entry 在 `_names_for` 中按 `scaffold_id` 运行时替换：

```python
# assembler.py `_names_for`（:290）
ring_stem = _ring_stem(numbered)          # 从 parent 的 stem_en/stem_zh 派生
if ring_stem:                             # 稠环/杂环 scaffold 词干
    entry = replace(entry, stem=ring_stem, coda="", omit_rule=..., aromatic=True)
elif sid == "carbocycle":
    entry = replace(entry, cyclic=True, ene_loc_omit=True, omit_rule=...)
sc_variant = (entry.variant or {}).get(sid)   # 苯环保留名等 scaffold 专属特例
```

### 4.2 特殊命名需求

如果你的环系需要特殊的命名处理，则可能需要：

- 在 `assembler.py` 的 `_names_for` 添加 worker（如 `_exocyclic_acid_names`/`_exocyclic_amide_names`/`_exocyclic_aldehyde_names`）
- 在 `assembler_prefixes.py` 中注册前缀构建规则
- 在 `chain_engine.py` 的 `_KIND_TABLE` entry 添加 `variant`（苯环保留名特例）

---

## 第五步: FG-环组合支持

FG-环组合的命名由以下路径承担：

- **苯系保留名**（benzoic/phenol/aniline/benzaldehyde/benzonitrile/benzamide/benzoate）：`chain_engine._KIND_TABLE` 各 entry 的 `variant["benzene"]`（要求 `scaffold_id=="benzene"` 且 `multiplicity==1`）
- **稠环/杂环 FG 收敛**：`express_ring_principal` 收敛为 FG 类别，词干由 `_ring_stem` 注入
- **环外酸**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_acid_names` → `cyclohexanecarboxylic acid`；多羧酸（multiplicity≥2）拼 …-di/tricarboxylic acid 且位次必带（P-65.2.2），苯单酸回落 benzoic acid 保留名
- **环外醛**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_aldehyde_names`（`assembler.py:170`）→ `cyclohexanecarbaldehyde`；多醛 → -dicarbaldehyde（P-66.6.1.1.3）；苯单醛仍走 benzaldehyde 保留名
- **环外酰胺**（`facts.relation == "exocyclic"`）：`assembler._exocyclic_amide_names` → `cyclohexanecarboxamide`
- **环二酸立体化学**：`layer4/cyclo_relative_stereo.py`（`acid` + `scaffold_id=="carbocycle"` → cis/trans）

---

## 文件改动清单

| Layer | 文件 | 改动内容 |
|-------|------|---------|
| L1 | `layer1/ring_systems.py` | 特殊环检测（通常无需改动） |
| L2 | `layer2/ring_scaffold.py` | 在 `_TEMPLATES` 新增 `{smiles, stem_en, stem_zh, naming_class}` |
| L2 | `layer2/kind_registry.py` | 通常无需改动（`_load_from_scaffold_specs` 自动注册） |
| L4 | `layer4/numbering_engine.py` | 通常无需改动（P-14.4 候选管线自动适用） |
| L5 | `layer5/assembler.py` | `_ring_stem` 词干注入 / exocyclic worker |
| L5 | `layer5/chain_engine.py` | `_KIND_TABLE` spec / `variant`（苯环保留名特例） |
| L5 | `layer5/assembler_prefixes.py` | 前缀构建规则（如需） |
| — | `tests/` | 添加 SMILES 测试用例 |

---

## 工作示例: 新增吡啶（pyridine）支持

以下是 NamePredict 中吡啶支持的实现路径（实际代码），理解它即可掌握新增保留名单杂环的完整流程。

### Step 1 — 环检测

吡啶的六元芳环（5C + 1N）被 SSSR 自动检测。未对 Layer1 做额外改动。

### Step 2 — Layer2 骨架识别

pyridine 走**数据驱动**识别，全部在 `ring_scaffold.py`：

1. **模板注册**：`_TEMPLATES` 已有 `pyridine` 条目（`stem_en="pyridine"`, `stem_zh="吡啶"`, `naming_class="monohetero"`），`_spec_from_template` 派生 ScaffoldSpec（retained=True）
2. **子图同构匹配**：`match_retained` 用 SMILES 模板做元素标注的子图同构，`resolve_ring_scaffold` 标注 `scaffold_id="pyridine"`

### Step 3 — 编号

无需 orienter。`numbering_engine` 的 P-14.4 候选管线对环枚举 2n 个候选；`_ring_hetero_start`（`numbering_engine.py`）把单杂环的唯一杂原子固定为 1 号位（P-14.4，pyridine 的 N=1）。

### Step 4 — 命名

`assembler._names_for` 查 `_KIND_TABLE`（`alkane` entry），按 `scaffold_id=="pyridine"` 从 `_ring_stem` 注入词干 `("pyridine","吡啶")`。FG 变体（pyridinecarboxylic 等）由 L2 `express_ring_principal` 收敛为 FG 类别后同样走 chain_engine。

### Step 5 — 测试

SMILES 测试用例示例：
- `c1ccncc1` — pyridine / 吡啶
- `c1ccncc1C(=O)O` — pyridine-2-carboxylic acid（吡啶甲酸 / 吡啶-2-甲酸）
- `c1ccncc1N` — pyridin-2-amine（吡啶-2-胺）

---

## 常见问题

### Q: 有 ScaffoldSpec 但没有骨架识别路径？

`ScaffoldSpec` 仅提供 stem 和编号策略的元数据；分子仍须在 `resolve_ring_scaffold` 的识别路径（`get_identity` 直查 / `match_retained` 子图同构）中被命中，其骨架才会使用该 scaffold。

### Q: 为什么有些环系没有独立的 .py 文件？

部分简单环系（如 furan, thiophene, pyrrole）不设独立模块——它们经 `ring_scaffold.py` 的 `_TEMPLATES` 批量注册，共享相同的编号策略，仅杂原子与 stem 不同。

### Q: 多环稠合系统如何处理？

多环稠合系统（如 naphthalene, indole）经 `ring_systems.py` 的 ring_system 分组识别，骨架身份由 `resolve_ring_scaffold` 解析（`match_retained` 子图同构）。编号方面，L2 注入 `numbering_scaffold` 事实（`{scaffold_id, labels}`，由 `numbering_scaffold_facts` 生成），L4 校验其存在性。

### Q: ScaffoldSpec 和 KindMeta 的关系？

`ScaffoldSpec` 是 stem 和编号的单一权威来源（Single Authority）。`kind_registry.py` 在 bootstrap 的唯一一步（`_load_from_scaffold_specs`，`:101`）读取所有 `ScaffoldSpec` 并转换为 `KindMeta` 注册。如果一个 kind 同时有 ScaffoldSpec 和手动 KindMeta 注册，ScaffoldSpec 覆盖前者（因为最后执行）。

> **源:** `src/namepredict/layer2/kind_registry.py:101`

---

## 相关文档

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器完整架构
- [[architecture/layer4-numbering]] — Layer4 编号引擎（P-14.4 候选管线）
- [[architecture/layer5-name-assembly]] — Layer5 名称组装
- [[concepts/functional-group-priority]] — FG 优先级与 IUPAC P-44 规则
