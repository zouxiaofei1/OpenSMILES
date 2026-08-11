# 新增环系指南

本文档分步说明如何为 NamePredict 流水线新增一个环系（ring system）——无论是新的 fused heterocycle（稠杂环）、bridged system（桥环）、spiro system（螺环），还是新的单杂环。

---

## 前置知识

阅读本指南前，建议先了解以下架构文档：

- [[architecture/overview]] — 六层流水线总览
- [[architecture/layer2-parent-selector]] — Layer2 母体选择器的核心架构（候选生成 + 评分排序）
- [[architecture/layer4-numbering]] — Layer4 定向与编号引擎
- [[architecture/layer5-name-assembly]] — Layer5 中英双语名称组装
- [[concepts/functional-group-priority]] — FG 优先级规则

---

## 流程概览

新增一个环系涉及五层流水线的改动，按顺序为：

1. **Layer1 环检测**（通常无需改动，SSSR 已自动识别）
2. **Layer2 母体选择**（主要工作量所在——添加 try 函数 + 注册 KindMeta / ScaffoldSpec）
3. **Layer4 编号定向**（添加 orienter + 注册到 `_kind_orienters()`）
4. **Layer5 名称组装**（添加命名函数 + 注册到 assembler dispatch chain）
5. **FG-环组合支持**（如环上有 COOH/CHO/CN/OH/NH2 等 FG 变体）

---

## 第一步: 环检测前提（Layer1）

大多数环系已被 Layer1 的 `ring_systems.py` 通过 **SSSR（Smallest Set of Smallest Rings）+ Union-Find 融合** 自动检测，无需额外工作。

> **源:** `E:\chem\src\namepredict\layer1\ring_systems.py:12-38`

SSSR 利用 RDKit 的 `GetRingInfo().AtomRings()` 获取所有最小环，然后通过 fusion graph（共享 >=2 原子的边）将稠合环合并为 ring system。桥环（共享 >=3 原子）和螺环（共享 1 原子）也被识别。

如果你的环系需要**超越 SSSR 的特殊识别**（例如识别特定的环取代模式、识别非 SSSR 的大环骨架），则需扩展 ring_systems.py 或创建新的检测模块。

### 环指纹（Ring Fingerprint）

`ring_fingerprint.py` 为每个环系生成指纹，格式为：

```
topology|sizes|fusion|hetero|aromatic
```

> **源:** `E:\chem\src\namepredict\layer1\ring_ir.py:29-34`

例如，吡啶（pyridine）的指纹类似 `mono|6||N|True`。Layer2 的 try 函数可以利用这些指纹快速筛选候选分子，而不必每次重新遍历原子。

---

## 第二步: 母体选择（Layer2）—— 主要工作量

这是整个流程中主要工作量所在。需要完成：注册 `ScaffoldSpec`（stem + 编号策略）、（可选）保留拓扑识别、以及特殊环系的薄层辅助。

### 2.1 骨架识别入口

母体选择统一走 P-44 规则管线（`rule_driven_parent_candidates`），环骨架身份由
`scaffold/ring_scaffold.resolve_ring_scaffold`（`ring_scaffold.py:37-45`）在表达阶段解析。
新环系的识别路径取决于其是否可拓扑匹配：

1. **ScaffoldSpec 直查** — 在 `scaffold/specs.py` 的规格表声明 `ScaffoldSpec`，`get_identity` 直接命中
2. **保留拓扑匹配** — 在 `scaffold/retained_registry.py` 的 `_TOPOLOGY` 添加拓扑条目（n_rings/n_atoms/hetero_Z/topology/aromatic），`match_systems` 遍历匹配
3. **薄层辅助模块** — 特殊环系（如 `scaffold/naphthalene.py`）在表达阶段补字段（`_naph_chains`/`_naph_parent_dict`）

参考现有实现：
- 单杂环数据驱动: `scaffold/specs.py` 的 `_monohetero`（`MONO_HETERO_SPECS`，`:174`）
- 保留拓扑: `scaffold/retained_registry.py` 的 `_TOPOLOGY`（benzene/pyridine/naphthalene/indole 4 条）
- 萘母体薄层: `scaffold/naphthalene.py`（桥环/螺环母体候选模块 `scaffold/polycyclic_parent.py` 已删除，L2 不产生 bridged/spiro 候选）

### 2.2 接入骨架识别（resolve_ring_scaffold）

`scaffold/ring_scaffold.py` 的 `resolve_ring_scaffold` 集中识别环骨架身份：

1. **specs 直查**：`scaffold/specs.py` 的 `get_identity(skeleton.scaffold_id)`（`:246`）直接命中 `ScaffoldSpec`
2. **保留拓扑匹配**：`_matched_id`（`ring_scaffold.py:10-13`）遍历 `retained_registry.match_systems(info)`（`:84`），找 `atom_ids` 与骨架原子集相等的系统
3. **兜底**：全碳环返回 `ScaffoldIdentity("carbocycle", ...)`（`_generic_carbocycle`）

> **源:** `E:\chem\src\namepredict\layer2\scaffold\ring_scaffold.py:37-45`

### 2.3 注册 ScaffoldSpec（新的命名骨架）

如果环系是 retained-name ring（如吡啶、萘、吲哚）或新的 scaffold，需在 `E:\chem\src\namepredict\layer2\scaffold\specs.py` 中添加一个 `ScaffoldSpec`。

`ScaffoldSpec` 字段：

| 字段 | 类型 | 说明 |
|------|------|------|
| `id` | `str` | 与 parent dict 的 `kind` 一致 |
| `naming_class` | `str` | 命名类别，如 `"monohetero"`, `"fused56"`, `"carbocycle"` |
| `stem_en` / `stem_zh` | `str` | 环系骨架的英文/中文 stem |
| `n_rings` | `int` | 环数（单环=1, 双环稠合=2, 三环=3） |
| `ring` | `str` | `"hetero"` 或 `"carbo"` |
| `retained` | `bool` | 是否为 IUPAC 保留名 |
| `fg_rank` | `int` | 官能团优先级（0 表示纯母体） |
| `numbering` | `NumberingPolicy` | 编号策略 |

`ScaffoldSpec` 中 `numbering` 的 `NumberingPolicy` 决定编号方式：

| mode | 适用场景 |
|------|---------|
| `"fixed_hetero"` | 单杂环，杂原子固定在 1 号位（吡啶、呋喃等） |
| `"fused56_fixed"` | 5+6 稠杂环，使用 1…7a 编号 |
| `"naph_family"` | 6+6 稠环，使用 1…8a 编号（萘、喹啉等） |
| `"carbocycle_free"` | 碳环自由编号（环烷烃） |
| `"anthracene_fixed"` | 蒽系三环编号（蒽、蒽醌） |
| `"fixed_roles"` | 固定角色编号（特定稠环） |

`standard_path` 指定固定的 locant 标签序列，例如 `FUSED56_LABELS = ("1","2","3","3a","4","5","6","7","7a")`。

> **源:** `E:\chem\src\namepredict\layer2\scaffold\specs.py:19-27, 31-38`

ScaffoldSpec 会自动被 `kind_registry.py` 的 bootstrap 加载为 `KindMeta`，提供 scoring 和 stem 元数据。

> **源:** `E:\chem\src\namepredict\layer2\kind_registry.py:207-214`

### 2.4 手动注册 KindMeta（无 ScaffoldSpec 的 FG 变体）

如果你的环系有 FG 组合变体（如 `pyridinol`, `pyridinamine`），但不是独立的 scaffold，则需要同时在 `kind_registry.py` 中手动注册 `KindMeta`。

方式一：加入 `_MISC_RING_FG` table（适合简单的 fg_rank 注册）：

> **源:** `E:\chem\src\namepredict\layer2\kind_registry.py:76-81`

方式二：在对应的 bootstrap 函数中调用 `_add()`：

```python
_add("pyridinol", fg=5, ret=True)
```

---

## 第三步: 编号定向（Layer4）

### 3.1 固定编号环系

对于有固定 IUPAC 编号的环系（稠环芳烃、杂环），添加 orienter 函数到 `E:\chem\src\namepredict\layer4\numbering.py`。

orienter 函数签名为：

```python
def _orient_{ring}(chain: list[int], parent: dict, substituents: list) -> list[int]:
```

它接收原始的原子顺序、parent dict、取代基列表，返回定向后的原子顺序。

常见的 orienter 模式：

- `_orient_pyridine` — 杂原子固定为 1 号位，然后选择使取代基位次最小的方向
  > **源:** `E:\chem\src\namepredict\layer4\numbering.py:199`
- `_orient_imidazole` — NH 固定为 1 号位，另一 N 取较小位次
  > **源:** `E:\chem\src\namepredict\layer4\numbering.py:205`
- `_orient_indole` — 稠环骨架已有固定顺序（从 ScaffoldSpec 的 standard_path），保持不变
  > **源:** `E:\chem\src\namepredict\layer4\numbering.py:231`

### 3.2 注册 orienter

将 orienter 注册到对应的字典中。系统通过以下调用链分发：

```
_kind_orienters() → _hetero_orienters() → _arene_orienters() → _aza_orienters()
```

对于氮杂环，在 `_aza_orienters()` 中注册：

> **源:** `E:\chem\src\namepredict\layer4\numbering.py:233`

对于稠环，在 `_fused_orienters()` 中注册：

> **源:** `E:\chem\src\namepredict\layer4\numbering.py:248`

对于其他环系，在 `_arene_orienters()` 或 `_hetero_orienters()` 中添加。

### 3.3 自由编号环系

对于无固定编号的环系（环烷烃、桥环、螺环），编号通过 constraint-based engine 计算：

- **简单碳环**: `_orient_cycloalkane`（取代基最低位次）
- **桥环**: `_orient_bridged`（von Baeyer 编号：桥头→最长桥→次长桥）
  > **源:** `E:\chem\src\namepredict\layer4\numbering.py:299`
- **复杂多约束**: 使用 `E:\chem\src\namepredict\layer4\locants\engine.py` 的 `choose_numbering`

### 3.4 FG locant 元组

如果环系有 OH/NH2 等 FG 变体且需要特殊 locant 提取逻辑，需更新 `_OH_KINDS` 和 `_AMINE_KINDS` 元组：

> **源:** `E:\chem\src\namepredict\layer4\numbering.py:412-413`

---

## 第四步: 名称组装（Layer5）

### 4.1 保留名环系

对于保留名环系（苯、吡啶、萘等），在 `E:\chem\src\namepredict\layer5\benzene_names.py` 中添加命名函数。

命名函数签名为：

```python
def {ring}_kind_names(kind: str, numbered: dict, build_prefix) -> tuple[str, str] | None:
```

- `kind` — 母体类型字符串
- `numbered` — Layer4 numbering 输出，包含 `parent`, `substituents`, locants 等
- `build_prefix` — 前缀构建函数（用于构建取代基前缀字符串）
- 返回 `(en_name, zh_name)` 或 `None`

对于吡啶系，`pyridine_kind_names` 作为分发中心处理 `pyridinecarboxylic`, `pyridinecarbonitrile`, `pyridinamine`, `pyridinol` 等子类：

> **源:** `E:\chem\src\namepredict\layer5\benzene_names.py:557`
> **源:** `E:\chem\src\namepredict\layer5\assembler.py:559`

### 4.2 系统命名环系

对于系统命名环系（环烷烃、桥环、螺环），命名逻辑通过 assembler 中的链式分发自动处理：

- 环烷烃 → `_cycloalkane_names` → `_ring_or_alkane`
- 桥环 → `_bridged_names` → `_ring_or_alkane`
- 螺环 → `_spiro_names` → `_ring_or_alkane`

> **源:** `E:\chem\src\namepredict\layer5\assembler.py:551-562`

大多数系统环系无需额外命名代码——stem 由 ScaffoldSpec 提供，assembler 自动拼接取代基前缀和后缀。

### 4.3 特殊命名需求

如果你的环系需要特殊的命名处理（如 fused ring locants "2a", "3a" 等），则可能需要：

- 在 `assembler_prefixes.py` 中注册前缀构建规则
- 在 `E:\chem\src\namepredict\layer5\stems.py` 中添加 stem 表

---

## 第五步: FG-环组合支持

如果环系需要支持 FG 变体（如吡啶甲酸 pyridinecarboxylic acid、吲哚甲醛 indolecarbaldehyde 等）：

### 5.1 FG 变体表达（principal 管线）

对于常见的 arene + FG 组合（naphthalenol、pyridinecarboxylic 等），由 layer2 的 principal 管线按 `kind` 组合处理（`principal_expression.py` 的 `_RETAINED_RING_KINDS`/`_resolved_ring_kind` + `kind_registry.py` 的 fg_rank 注册）。

### 5.2 注册

对于带 FG 变体的保留母环，在 `scaffold/specs.py` 登记对应 `ScaffoldSpec`（如 `pyrazolamine`/`thiazolamine` 走 `_monohetero`），或在 `kind_registry.py` 的 `_MISC_RING_FG` 注册环 FG 变体的 fg_rank（如 `cycloalkane_polycarboxylic`/`cycloalkanedione`），确保 FG-ring 组合的 `kind` 有对应 `KindMeta`，以获得正确的 `fg_rank` 评分。

---

## 文件改动清单

| Layer | 文件 | 改动内容 |
|-------|------|---------|
| L1 | `layer1/ring_systems.py` | 特殊环检测（通常无需改动） |
| L2 | `layer2/scaffold/specs.py` | 新增 `ScaffoldSpec`（stem + 编号策略） |
| L2 | `layer2/scaffold/retained_registry.py` | 新增保留拓扑条目（如可拓扑匹配） |
| L2 | `layer2/scaffold/{ring_name}.py` | 特殊环系薄层辅助模块（如需要） |
| L2 | `layer2/kind_registry.py` | 注册 `KindMeta`（如 FG 变体无 ScaffoldSpec） |
| L4 | `layer4/numbering.py` | 添加 orienter + 注册到 `_kind_orienters()` 链 |
| L4 | `layer4/locants/engine.py` | 约束模式（如需要复杂约束编号） |
| L5 | `layer5/benzene_names.py` 或 `{ring}_names.py` | 命名函数 + 分发注册 |
| L5 | `layer5/assembler.py` | 注册到 dispatch chain（如 `_ring_or_alkane`） |
| L5 | `layer5/assembler_prefixes.py` | 前缀构建规则（如需） |
| — | `tests/` | 添加 SMILES 测试用例 |

---

## 工作示例: 新增吡啶（pyridine）支持

以下是 NamePredict 中吡啶支持的实现路径（简化版，实际代码更完整）。以它为参考，可以理解新增一个保留名单杂环的完整流程。

### Step 1 — 环检测

吡啶的六元芳环（5C + 1N）被 SSSR 自动检测。未对 Layer1 做额外改动。

### Step 2 — Layer2 骨架识别

pyridine 走**数据驱动**识别：

1. **ScaffoldSpec**：`scaffold/specs.py` 的 `MONO_HETERO_SPECS` 已有
   `_monohetero("pyridine", "pyridine", "吡啶")`（`:175`，提供 stem + `fixed_hetero` 编号策略）
2. **保留拓扑匹配**：`scaffold/retained_registry.py` 的 `_TOPOLOGY` 已有 pyridine 条目
   （`n_rings=1, n_atoms=6, hetero_Z=(7,), topology="mono", aromatic=True`），
   `resolve_ring_scaffold` 经 `match_systems` 标注 `scaffold_id="pyridine"`
3. **FG 变体**：pyridinecarboxylic / pyridinamine / pyridinol 等 kind 由 `kind_registry`
   注册 fg_rank，命名走 L5 的 `pyridine_kind_names`（`benzene_names.py:315`）分发

特殊环系（非拓扑匹配）可新建薄层模块在表达阶段补字段（参考 `scaffold/naphthalene.py`）。

### Step 2b — ScaffoldSpec 声明

文件 `E:\chem\src\namepredict\layer2\scaffold\specs.py`:

```python
_monohetero("pyridine", "pyridine", "吡啶"),
```

这会生成 `ScaffoldSpec(id="pyridine", naming_class="monohetero", stem_en="pyridine", stem_zh="吡啶", n_rings=1, ring="hetero", retained=True, numbering=NumberingPolicy(mode="fixed_hetero"))`。

> **源:** `E:\chem\src\namepredict\layer2\scaffold\specs.py:175`

### Step 3 — 编号 orienter

文件 `E:\chem\src\namepredict\layer4\numbering.py`:

```python
def _orient_pyridine(chain, parent, substituents):
    return _orient_ring_fixed(chain, parent, substituents, "n_idx")
    # 将 N 原子固定在 1 号位, 然后选取代基位次最小的方向
```

注册到 `_aza_orienters()`:

```python
"pyridine": _orient_pyridine,
"pyridinecarboxylic": _orient_pyridinecarboxylic,
"pyridinamine": _orient_pyridinamine,
"pyridinol": _orient_pyridinol,
```

> **源:** `E:\chem\src\namepredict\layer4\numbering.py:199, 233`

### Step 4 — 命名

文件 `E:\chem\src\namepredict\layer5\benzene_names.py` 中的 `pyridine_kind_names`:

```python
def pyridine_kind_names(kind, numbered, build_prefix):
    # 根据 kind 分派到具体的命名函数:
    # pyridinecarboxylic → "pyridine-N-carboxylic acid" / "吡啶-N-甲酸"
    # pyridinamine → "pyridin-N-amine" / "吡啶-N-胺"
    # ...
```

文件 `E:\chem\src\namepredict\layer5\assembler.py` 中的 dispatch:

```python
top = pyridine_kind_names(kind, numbered, _build_prefix)
```

> **源:** `E:\chem\src\namepredict\layer5\assembler.py:559`

### Step 5 — 测试

SMILES 测试用例示例：
- `c1ccncc1` — pyridine / 吡啶
- `c1ccncc1C(=O)O` — pyridine-2-carboxylic acid（吡啶甲酸 / 吡啶-2-甲酸）
- `c1ccncc1N` — pyridin-2-amine（吡啶-2-胺）

---

## 常见问题

### Q: 有 ScaffoldSpec 但没有骨架识别路径？

`ScaffoldSpec` 仅提供 stem 和编号策略的元数据；分子仍须在 `resolve_ring_scaffold` 的识别路径（specs 直查 / 保留拓扑匹配）中被命中，其骨架才会使用该 scaffold。

### Q: 为什么有些环系没有独立的 .py 文件？

部分简单环系（如 furan, thiophene, pyrrole）不设独立模块——它们经 `scaffold/specs.py` 的 `_monohetero` 批量注册（`MONO_HETERO_SPECS`，`:174`），共享相同的 `fixed_hetero` 编号策略，仅杂原子与 stem 不同。

> **源:** `E:\chem\src\namepredict\layer2\scaffold\specs.py:174-198`

### Q: 多环稠合系统如何处理？

多环稠合系统（如 anthracene, naphthalene）经 `ring_systems.py` 的 ring_system 分组识别，骨架身份由 `resolve_ring_scaffold` 解析（保留拓扑匹配，或 `_poly_carbo`/`_naph`/`_fused56` ScaffoldSpec）。编号方面，`ScaffoldSpec` 的 `standard_path` 提供固定的 locant 标签（如 `FUSED56_LABELS` 或 `NAPH_LABELS`），Layer4 依据 `NumberingPolicy` 定向，无需额外 orienter。

### Q: ScaffoldSpec 和 KindMeta 的关系？

`ScaffoldSpec` 是 stem 和编号的单一权威来源（Single Authority）。`kind_registry.py` 在 bootstrap 的最后一步（`_load_from_scaffold_specs`，`:207`）读取所有 `ScaffoldSpec` 并转换为 `KindMeta` 注册。如果一个 kind 同时有 ScaffoldSpec 和手动 KindMeta 注册，ScaffoldSpec 覆盖前者（因为 `_load_from_scaffold_specs` 最后执行）。

> **源:** `E:\chem\src\namepredict\layer2\kind_registry.py:207-214, 224`

---

## 相关文档

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器完整架构
- [[architecture/layer4-numbering]] — Layer4 定向与编号引擎
- [[architecture/layer5-name-assembly]] — Layer5 名称组装
- [[concepts/functional-group-priority]] — FG 优先级与 IUPAC P-44 规则
