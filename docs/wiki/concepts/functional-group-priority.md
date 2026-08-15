# 官能团优先层级 (Functional Group Priority Hierarchy)

> **核心概念** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer5-name-assembly]], [[concepts/atom-ownership]]

---

## 概述

IUPAC 有机命名中，一个分子可能同时含有多种官能团（Functional Group, FG），例如羟基酸（含 -COOH 和 -OH）、氨基酮（含 -NH&#8322; 和 >C=O）、氰基酯（含 -CN 和 -COOR）等。按照 IUPAC P-41 规则，这些官能团之间存在严格的**优先顺序**：优先级最高的官能团成为 **principal characteristic group（主特征基团）**，以母体后缀（suffix）表达；较低优先级的官能团退化为取代基前缀（prefix）。

NamePredict 将这一优先级体系内建于 **`kind_registry.py` 的 `fg_rank` 字段**，并在三层流水线中接力使用：

1. **Layer 1（`analyzer.py`）**: 检测 FG 时采用排他性优先级（如 carboxyl 排斥 ester/amide/anhydride）
2. **Layer 2（`scoring.py`）**: 将 `fg_rank` 作为 P-44 评分 tuple 的第二维，决定母体选择
3. **Layer 5（`chain_engine.py`）**: 直接查 `_KIND_TABLE` 分派后缀（kind 收敛在 L2 `_chain_kind`，无 `typed_kinds` 模块），高优先级 FG 获得后缀，低优先级 FG 转为前缀

> 13 个扩展 FG（sulfoxide/sulfone/sulfonate/sulfonamide/sulfonic_acid/sulfonyl_chloride/phosphate/boronic/carbamate/carbonate/urea/guanidine/hydrazine）在 layer1 检测中不存在，随之退出 fg_rank 体系。

---

## fg_rank 优先级表

以下表格展示 NamePredict 中实现的官能团优先级，遵循 IUPAC P-41 顺序（数值越大优先级越高）。表格依据 `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY`（`compatibility_rank`）生成，经 `kind_registry._KIND_CLASS` 映射到 kind。

| `fg_rank` | 官能团种类 (kind) | 英文后缀示例 | 中文后缀示例 | 说明 |
|:---:|---|---|---|---|
| **14** | `acid` | -oic acid | -酸 | P-65.1: 羧酸为最高优先级链状 FG |
| **12** | `anhydride` | -oic anhydride | -酸酐 | 羧酸酐（registry 保留；`_CHAIN_FG` 无此类，链式 acid 无酸酐表达） |
| **11** | `ester` | -oate | -酸酯 | 羧酸衍生物（酯类） |
| **10** | `acyl_chloride`, `acyl_bromide` | -oyl chloride / -oyl bromide | -酰氯 / -酰溴 | 酰卤 |
| **9** | `amide` | -amide | -酰胺 | 酰胺 |
| **8** | `nitrile`, `isocyanate`, `isothiocyanate` | -nitrile / isocyanate | -腈 / 异氰酸酯 | C&#8801;N 和累积双键系统 |
| **7** | `aldehyde` | -al / carbaldehyde | -醛 / 甲醛 | 醛基 |
| **6** | `ketone` | -one | -酮 | 酮（`dione` 已不再由 L2 产生，二酮由 L5 chain_engine `mult_ok` 生成式命名） |
| **5** | `alcohol` | -ol | -醇 | 羟基 |
| **4** | `thiol` | -thiol | -硫醇 | 巯基 |
| **3** | `amine` | -amine | -胺 | 氨基 |
| **2** | `sulfide` | sulfide | 硫醚 | 低优先级杂原子 FG（LEGACY_COMPAT） |
| **0** | `ether`, `alkane`, `alkene`, `alkyne` | ether / -ane / -ene / -yne | 醚 / -烷 / -烯 / -炔 | 无 principal FG 的母体（醚始终为前缀） |

> **源:** `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY` (P-41 链状 FG), `src/namepredict/layer2/kind_registry.py` 的 `_principal_rank` 投影

> 注：数量派生 kind（`diacid`/`polycarboxylic`/`diol`/`triol`/`diamine`/`triamine`/`tetraamine`）与组合 kind（`cycloalcohol`/`cycloketone`/`cycloamine`/`cycloalkane_polycarboxylic` 等）已全部删除——链式 acid/alcohol/amine/ketone 对任意主基团数 kind 恒为基团名，数量由 `principal_expression_facts.multiplicity` 承载。苯/杂环保留名（benzoic/phenol/aniline 等）由 L5 chain_engine `variant` 提供，不占独立 kind 行。

---

## KindMeta 数据结构

`kind_registry.py:34-42` 定义了 `KindMeta` 数据类，作为所有母体种类的统一元数据容器：

```python
@dataclass(frozen=True)
class KindMeta:
    kind: str              # 母体种类标识符（如 "acid", "ketone", "benzene"）
    en: str | None         # 英文 stem 名称
    zh: str | None         # 中文 stem 名称
    fg_rank: int           # FG 优先级 (0-13)，0 表示无 principal FG
    ring: str              # 环系类型: "none" | "hetero" | "carbo"
    n_rings: int           # 环的数量
    retained: bool         # 是否为 IUPAC 保留名（retained name），影响 L2 评分
```

所有 scaffold 母体种类通过模块级 `_bootstrap()` 函数统一注册（`kind_registry.py:133-135`），注册**只有一步**：`_load_from_scaffold_specs()`（从 `ring_scaffold.all_specs()` 读取，`:123`）。**`ring_scaffold.py` 的 `_TEMPLATES`** 是 **stem 的最终权威来源**。链式 FG kind 不预先注册——无 `_KIND_CLASS`/`_load_chain_fg`/`all_kinds`，`fg_rank` 由 `principal.legacy_rank` 实时投影（`_principal_rank`，`kind_registry.py:18`）。

---

## 三层 FG 生命周期

### Layer 1: 检测（Detect）—— 排他性优先级

`src/namepredict/layer1/analyzer.py` 在检测 FG 时，必须处理**化学上的包含关系**：一个 carbon 如果已经是 carboxyl carbon，就不能同时被识别为 ester 或 amide。为实现这一点，Layer 1 采用"排他性检测"模式：

- **Carboxyl 优先于 ester/amide/anhydride**: `_is_carboxyl_carbon()` 匹配 C(=O)OH 或 C(=O)O⁻ 模式，ester/amide 检测显式判定酸性氧邻居为 False 后才匹配。
- **Anhydride 桥氧从 ether 中排除**: `_is_ether_oxygen()` 先检查 `_is_anhydride_bridge_o()`。

13 个扩展 FG（carbamate/carbonate/urea/guanidine/sulfoxide/sulfone/...）在 layer1 中不存在，因此无跨模块排他链（carbamate 排除 ester、urea 排除 amide 等）——排他检测只发生在 analyzer.py 内部的核心羰基族。共享的羰基检测原语集中在 `_carbonyl_common.py`。

所有检测到的 FG 被汇总为结构化列表（20 个 FG 列表键 + 18 个布尔标志）和类型化的 `fg_inventory`，构成 info dict 传递给 Layer 2。

> **源:** `src/namepredict/layer1/analyzer.py`, `src/namepredict/layer1/_carbonyl_common.py` | 详情见 [[architecture/layer1-analyzer]]

### Layer 2: 评分与选择（Score & Select）

`src/namepredict/layer2/scoring.py:49-58` 定义了 11 维 P-44 评分 tuple，`fg_rank` 经 `principal_group_class` 编码在第二维：

```python
(principal_group_class,   # FG 类别 rank（principal contract，rank=0 即无主官能团）
 principal_group_count,   # 主官能团实例数
 sides_ok,                # 0/1 — 侧链是否全部表达
 is_hetero_ring,          # 0/1 — 杂环加分
 is_carbo_ring,           # 0/1 — 碳环加分
 n_rings,                 # 环数
 ring_size,               # 环原子数
 retained_bonus,          # 保留名加分
 n_unsat,                 # 不饱和键数（字段驱动判读）
 n_carbons,               # 碳原子数
 -n_unhandled_side)       # 未识别侧链惩罚
```

这一评分体系使得含羧酸的链状母体（`fg_rank=14`）必然优先于含酮的母体（`fg_rank=6`），无论链长或环数如何。主官能团有无由 `principal_group_class` 承载（rank=0 即无主官能团），无独立 `has_principal_fg` 布尔位。

母体候选的生成以 **P-44 规则驱动管线**为主（`rule_driven_parent_candidates`，见
[[architecture/layer2-parent-selector]]）。`fg_rank` 经 `kind_registry._KIND_CLASS` 映射到 FG 枚举后由 `principal.legacy_rank` 投影（P-41 `compatibility_rank` 为单一权威）。

> **源:** `src/namepredict/layer2/scoring.py:49-58`, `kind_registry.py:11-24` | 详情见 [[architecture/layer2-parent-selector]]

### Layer 5: 后缀分派（Suffix Dispatch）

Layer5 由 `_names_for`（`assembler.py:101`）查 **`chain_engine._KIND_TABLE`**（10 个 `_Chain` spec）渲染词干、不饱和段、位次与环前缀（kind 收敛在 L2 `_chain_kind`，无 `typed_kinds` 模块）。`_Chain.variant` 字段按 multiplicity 切换数量后缀（alcohol→diol/triol/tetraol，amine→diamine/triamine/tetraamine，acid→dioic acid）。特殊 case 走 worker（`_exocyclic_acid_names`/`_exocyclic_amide_names`/`_parent_stem_names`）。

调度是分层的：L2 收敛 kind → L5 查链引擎 → 命中后直接返回。这确保了被选为母体的 principal FG 获得后缀，而劣后 FG 在 Layer 3 中被转为取代基前缀（如 hydroxy-、oxo-、amino-）。

> **源:** `src/namepredict/layer5/chain_engine.py`, `src/namepredict/layer5/assembler.py:101-136` | 详情见 [[architecture/layer5-name-assembly]]

---

## 互斥排除（P-41 互斥约束）

IUPAC P-41 规定某些 FG 之间不能作为母体共存。NamePredict 通过 **`principal.py` 的 `select_principal_group`** 结构性实现互斥：P-44 只选单个最高优先级主官能团（`min(eligible, key=priority)`），低优先级 FG 一律成为取代基，不需要逐候选互斥检查。

> 无 `fg_helpers.py` 与 `_no_fgs(info, keys)` 互斥谓词——互斥由 `select_principal_group` 结构性单选择实现。

---

## 注册体系：kind_registry 作为元数据聚合层

`src/namepredict/layer2/kind_registry.py` 是 FG 元数据的**注册与查询中心**：链状 FG 的 `fg_rank` 经 `principal.legacy_rank` 从 `PRINCIPAL_REGISTRY.compatibility_rank` 实时投影（P-41 单一权威，`_principal_rank`），保留 scaffold 的 stem 由 `ring_scaffold.py` 的 `_TEMPLATES` 派生（bootstrap 唯一一步 `_load_from_scaffold_specs`）。所有对 FG 优先级的查询通过此模块的公共 API 进行：

| API | 功能 |
|---|---|
| `kind_registry.fg_rank(kind)` | 返回 kind 的 `fg_rank` 值 |
| `kind_registry.has_principal_fg(kind)` | 是否有 principal FG（`fg_rank > 0`） |
| `kind_registry.is_hetero_ring(kind)` / `is_carbo_ring(kind)` | 环系类型查询 |
| `kind_registry.retained_bonus(kind)` | 是否为保留名 |
| `kind_registry.parent_names(kind)` | 返回 `(en_stem, zh_stem)` 元组 |

新增 FG 种类时，若走 principal typed 管线（acid/alcohol/ketone/amine 等），在 `principal.py` 的 `PRINCIPAL_REGISTRY` 注册 `PrincipalFeatureSpec`（SUFFIX 档），并在 `principal_expression.py` 的 `_CHAIN_KINDS` 表达表中声明 kind 即可。评分、候选收集与组装逻辑不变，仍保持**开放-封闭原则**。

> **源:** `src/namepredict/layer2/kind_registry.py`

---

## 总结

官能团优先级是 NamePredict 命名正确性的基石。从 Layer 1 的排他性检测、Layer 2 的 `fg_rank` 评分、Layer 5 的后缀分派，整个系统严格遵循 IUPAC P-41 和 P-44 的层级化思维：

1. **检测层保证"一碳一 FG"**——排他性模式使得同一碳原子不会被重复识别为多个 FG
2. **评分层保证"优先级强制"**——`fg_rank` 在 P-44 tuple 中高于环数和碳链长度，确保羧酸绝对优先于酮
3. **组装层保证"劣后 FG 降级"**——低优先级 FG 从候选后缀降为取代基前缀
4. **互斥层保证"互斥正确"**——`select_principal_group` 的单选择结构性排除无法共存的 FG 组合
