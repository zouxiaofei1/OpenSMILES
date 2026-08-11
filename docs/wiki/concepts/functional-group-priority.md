# 官能团优先层级 (Functional Group Priority Hierarchy)

> **核心概念** | 关联: [[architecture/overview]], [[architecture/layer1-analyzer]], [[architecture/layer2-parent-selector]], [[architecture/layer5-name-assembly]], [[concepts/atom-ownership]]

---

## 概述

IUPAC 有机命名中，一个分子可能同时含有多种官能团（Functional Group, FG），例如羟基酸（含 -COOH 和 -OH）、氨基酮（含 -NH&#8322; 和 >C=O）、氰基酯（含 -CN 和 -COOR）等。按照 IUPAC P-41 规则，这些官能团之间存在严格的**优先顺序**：优先级最高的官能团成为 **principal characteristic group（主特征基团）**，以母体后缀（suffix）表达；较低优先级的官能团退化为取代基前缀（prefix）。

NamePredict 将这一优先级体系内建于 **`kind_registry.py` 的 `fg_rank` 字段**，并在三层流水线中接力使用：

1. **Layer 1（`analyzer.py`）**: 检测 FG 时采用排他性优先级（如 carboxyl 排斥 ester/amide/anhydride）
2. **Layer 2（`scoring.py`）**: 将 `fg_rank` 作为 P-44 评分 tuple 的第二维，决定母体选择
3. **Layer 5（`assembler.py`）**: 根据 parent kind 分派后缀，高优先级 FG 获得后缀，低优先级 FG 转为前缀

---

## fg_rank 优先级表

以下表格展示 NamePredict 中实现的官能团优先级，遵循 IUPAC P-41 顺序（数值越大优先级越高）。表格依据 `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY`（`compatibility_rank`）生成，经 `kind_registry._KIND_CLASS` 映射到 kind，同时合并了从 `_ARENE_NAMED`、`_MISC_RING_FG` 等注册表注入的环系 FG 种类。

| `fg_rank` | 官能团种类 (kind) | 英文后缀示例 | 中文后缀示例 | 说明 |
|:---:|---|---|---|---|
| **14** | `acid`, `diacid`, `polycarboxylic` | -oic acid / -dioic acid | -酸 / -二酸 | P-65.1: 羧酸为最高优先级链状 FG |
| **14** | `benzoic`, `furancarboxylic`, `pyridinecarboxylic` 等 | benzoic acid / ... | 苯甲酸 / ... | 芳环/杂环羧酸保留名，同属 carboxylic acid 类 |
| **13** | `sulfonic_acid` | sulfonic acid | 磺酸 | S 酸类，略低于羧酸 |
| **12** | `anhydride`, `boronic` | -oic anhydride / boronic acid | -酸酐 / 硼酸 | 羧酸酐与 B 酸 |
| **11** | `ester`, `diester`, `carbamate`, `carbonate` | -oate / -dioate | -酸酯 / 氨基甲酸酯 | 羧酸衍生物（酯类） |
| **11** | `sulfonate` | sulfonate | 磺酸酯 | 磺酸酯 |
| **10** | `acyl_chloride`, `acyl_bromide` | -oyl chloride / -oyl bromide | -酰氯 / -酰溴 | 酰卤 |
| **10** | `sulfonyl_chloride` | sulfonyl chloride | 磺酰氯 | 磺酰卤 |
| **9** | `amide`, `urea`, `guanidine`, `sulfonamide` | -amide / urea | -酰胺 / 脲 | 酰胺和 N-衍生物 |
| **8** | `nitrile`, `isocyanate`, `isothiocyanate` | -nitrile / isocyanate | -腈 / 异氰酸酯 | C&#8801;N 和累积双键系统 |
| **7** | `aldehyde` | -al / carbaldehyde | -醛 / 甲醛 | 醛基 |
| **6** | `ketone`, `dione`, `cycloketone`, `sulfone` | -one / -dione | -酮 / -二酮 | 酮（含砜，S-氧化态等同酮优先级） |
| **5** | `alcohol`, `diol`, `triol`, `cycloalcohol` | -ol / -diol / -triol | -醇 / -二醇 / -三醇 | 羟基 |
| **4** | `thiol`, `hydrazine` | -thiol | -硫醇 | 巯基 / 肼 |
| **3** | `amine`, `diamine`, `triamine`, `tetraamine`, `cycloamine`, `aniline` | -amine / -diamine | -胺 / -二胺 | 氨基 |
| **3** | `sec_amine`, `tert_amine` | N-...-amine | N-...-胺 | 二级/三级胺 |
| **2** | `phosphate`, `phosphonic` | phosphate / phosphonic acid | 磷酸酯 / 膦酸 | P-官能团 |
| **2** | `sulfide`, `sulfoxide` | sulfide / sulfoxide | 硫醚 / 亚砜 | 低优先级杂原子 FG |
| **0** | `ether`, `alkane`, `alkene`, `alkyne`, `cycloalkane`, `bridged` | ether / -ane / -ene / -yne | 醚 / -烷 / -烯 / -炔 | 无 principal FG 的母体（醚始终为前缀） |

> **源:** `src/namepredict/layer2/principal.py` 的 `PRINCIPAL_REGISTRY` (P-41 链状 FG), `src/namepredict/layer2/kind_registry.py` 的 `_ARENE_NAMED` (芳烃命名 FG), `_MISC_RING_FG` (杂项环 FG)

---

## KindMeta 数据结构

`kind_registry.py:45-53` 定义了 `KindMeta` 数据类，作为所有母体种类的统一元数据容器：

```python
@dataclass(frozen=True)
class KindMeta:
    kind: str              # 母体种类标识符（如 "acid", "ketone", "pyridinecarboxylic"）
    en: str | None         # 英文 stem 名称
    zh: str | None         # 中文 stem 名称
    fg_rank: int           # FG 优先级 (0-13)，0 表示无 principal FG
    ring: str              # 环系类型: "none" | "hetero" | "carbo"
    n_rings: int           # 环的数量
    retained: bool         # 是否为 IUPAC 保留名（retained name），影响 L2 评分
```

所有母体种类通过模块级 `_bootstrap()` 函数统一注册（`kind_registry.py:217-227`），注册顺序为链状 FG → 芳烃 FG → 环 FG 变体 → 环烷 → bridged → 饱和杂环 → ScaffoldSpec。最后加载的 `_load_from_scaffold_specs()` 来自 `scaffold/specs.py`，是 **stem 的最终权威来源**——如果已有注册条目，ScaffoldSpec 会覆盖。

---

## 三层 FG 生命周期

### Layer 1: 检测（Detect）—— 排他性优先级

`src/namepredict/layer1/analyzer.py` 在检测 FG 时，必须处理**化学上的包含关系**：一个 carbon 如果已经是 carboxyl carbon，就不能同时被识别为 ester 或 amide。为实现这一点，Layer 1 采用"排他性检测"模式：

- **Carboxyl 优先于 ester/amide/anhydride**: `_is_carboxyl_carbon()` (line 60) 匹配 C(=O)OH 或 C(=O)O⁻ 模式。`_is_ester_carbon()` (line 194) 显式判定 `_has_acid_o_neighbor(atom)` 为 False 后才匹配；`_is_amide_carbon()` (line 89) 同样排斥 `_has_acid_o_neighbor(atom)`。
- **Carbamate 优先于 ester**: `_is_ester_carbon()` 调用 `_is_carbamate_carbon()`，若碳原子属于 carbamate（O=C-O-N 模式），则排除在 ester 之外。
- **Urea 优先于 amide**: 同理，`_is_amide_carbon()` 调用 `_is_urea_carbon()` 排除 urea 碳。
- **Anhydride 桥氧从 ether 中排除**: `_is_ether_oxygen()` (line 127) 先检查 `_is_anhydride_bridge_o()`。

所有检测到的 FG 被汇总为结构化列表（如 `carboxyls: [{"c_idx": 5, "anion": False}, ...]`）和布尔标志（`has_acid: True`），构成 info dict 传递给 Layer 2。

> **源:** `src/namepredict/layer1/analyzer.py:60-63, 89-95, 194-201` | 详情见 [[architecture/layer1-analyzer]]

### Layer 2: 评分与选择（Score & Select）

`src/namepredict/layer2/scoring.py:73-81` 定义了 11 维 P-44 评分 tuple，其中 `fg_rank` 是第二维（仅次于 `has_principal_fg` 布尔标志）：

```python
(has_principal_fg,    # 0/1 — 是否有主特征基团
 fg_class_rank,       # ← fg_rank 值 (0-13)，数值越大越优先
 sides_ok,            # 0/1 — 侧链是否全部表达
 is_hetero_ring,      # 0/1 — 杂环加分
 is_carbo_ring,       # 0/1 — 碳环加分
 n_rings,             # 环数
 ring_size,           # 环原子数
 retained_bonus,      # 保留名加分
 n_unsat,             # 不饱和键数
 n_carbons,           # 碳原子数
 -n_unhandled_side)   # 未识别侧链惩罚
```

这一评分体系使得含羧酸的链状母体（`fg_rank=13`, `has_principal_fg=1`）必然优先于含酮的母体（`fg_rank=6`），无论链长或环数如何。

母体候选的生成以 **P-44 规则驱动管线**为主（`rule_driven_parent_candidates`，见
[[architecture/layer2-parent-selector]]），评分由 `scoring.py` 的 11 维 tuple 仲裁。`fg_rank` 经
`kind_registry._KIND_CLASS` 映射到 FG 枚举后由 `principal.legacy_rank` 投影（P-41 `compatibility_rank`
为单一权威），不存在任何 try 函数注册表。

> **源:** `src/namepredict/layer2/scoring.py:65-68`, `kind_registry.py:10-33` | 详情见 [[architecture/layer2-parent-selector]]

### Layer 5: 后缀分派（Suffix Dispatch）

`src/namepredict/layer5/assembler.py:493-499` 的 `_names_for()` 函数是名称组装的中枢调度器。它根据 parent kind 分层分发：

1. **`special_fg_names()`**: 处理特殊 FG 类别（酯的 alkoxy 配对、sulfonate、carbamate 等非标准后缀）
2. **`p_fg_names()`**: 处理 phosphate/phosphonic 的 P-FG 后缀
3. **`_hetero_names()`**: 分发 alcohol/thiol/ether/sulfide/amine 类后缀
4. **`_carbonyl_names()`**: 分发 acid/aldehyde/amide/nitrile/ester/ketone 类后缀
5. **`_unsat_or_alkane()`**: 无 principal FG 时的烃类后缀（-ane/-ene/-yne）

调度逻辑是分层的：先尝试最高优先级路径，命中后直接返回。这确保了被选为母体的 principal FG 获得后缀，而劣后 FG 在 Layer 3 中被转为取代基前缀（如 hydroxy-、oxo-、amino-）。

> **源:** `src/namepredict/layer5/assembler.py:493-499`, `src/namepredict/layer5/stems.py:1-376` | 详情见 [[architecture/layer5-name-assembly]]

---

## 互斥排除（P-41 互斥约束）

IUPAC P-41 规定某些 FG 之间不能作为母体共存——当某个候选母体被选中时，分子中不得存在与其冲突的更高或同级 FG。NamePredict 通过**互斥检查谓词 `_no_fgs(info, keys)`**（`fg_helpers.py`）实现：`keys` 是一组 `has_*` 布尔键，若其中任何一个为 `True`，则该候选母体被排除。

各候选的互斥 keys 由调用方**内联**传入 `_no_fgs(info, keys)`。`_no_fgs`（`fg_helpers.py:26`）是唯一的共享互斥原语；每类母体的排斥集合内联在各自模块的 `keys` 参数里。这是一种声明式的约束表达：**`keys` 元组定义了"该母体不容忍的 FG 集合"**。

> **源:** `src/namepredict/layer2/fg_helpers.py:26`

---

## 注册体系：kind_registry 作为元数据聚合层

`src/namepredict/layer2/kind_registry.py` 是 FG 元数据的**注册与查询中心**：链状 FG 的 `fg_rank` 经
`_KIND_CLASS` 映射到 FG 枚举后由 `principal.PRINCIPAL_REGISTRY` 的 `compatibility_rank` 投影（P-41
单一权威），保留 scaffold 的 stem 由 `scaffold/specs.py` 提供（bootstrap 最后加载覆盖）。所有对 FG
优先级的查询通过此模块的公共 API 进行：

| API | 功能 |
|---|---|
| `kind_registry.fg_rank(kind)` | 返回 kind 的 `fg_rank` 值 |
| `kind_registry.has_principal_fg(kind)` | 是否有 principal FG（`fg_rank > 0`） |
| `kind_registry.is_hetero_ring(kind)` / `is_carbo_ring(kind)` | 环系类型查询 |
| `kind_registry.retained_bonus(kind)` | 是否为保留名 |
| `kind_registry.parent_names(kind)` | 返回 `(en_stem, zh_stem)` 元组 |

新增 FG 种类时，若走 principal typed 管线（acid/alcohol/ketone/amine 等），在 `principal.py` 的
`PRINCIPAL_REGISTRY` 注册 `PrincipalFeatureSpec`（SUFFIX 档），并在 `principal_expression.py`
的 `_RETAINED_RING_KINDS` / `_CHAIN_KINDS` 等表达表中声明 kind 与字段即可。评分、候选收集与组装
逻辑不变，仍保持**开放-封闭原则**。

> **源:** `src/namepredict/layer2/kind_registry.py:89-227`

---

## 总结

官能团优先级是 NamePredict 命名正确性的基石。从 Layer 1 的排他性检测、Layer 2 的 `fg_rank` 评分、Layer 5 的后缀分派，到互斥检查规则，整个系统严格遵循 IUPAC P-41 和 P-44 的层级化思维：

1. **检测层保证"一碳一 FG"**——排他性模式使得同一碳原子不会被重复识别为多个 FG
2. **评分层保证"优先级强制"**——`fg_rank` 在 P-44 tuple 中高于环数和碳链长度，确保羧酸绝对优先于酮
3. **组装层保证"劣后 FG 降级"**——低优先级 FG 从候选后缀降为取代基前缀
4. **互斥层保证"互斥正确"**——无法共存的 FG 组合通过 `_no_fgs` 互斥检查（`fg_helpers.py`）在候选阶段即被排除
