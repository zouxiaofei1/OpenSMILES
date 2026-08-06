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

以下表格展示 NamePredict 中实现的官能团优先级，遵循 IUPAC P-41 顺序（数值越大优先级越高）。表格依据 `src/namepredict/layer2/kind_registry.py:21-36` 中的 `_CHAIN_FG` 注册表生成，同时合并了从 `_ARENE_NAMED`、`_H5_COOH`、`_MISC_RING_FG` 等注册表注入的环系 FG 种类。

| `fg_rank` | 官能团种类 (kind) | 英文后缀示例 | 中文后缀示例 | 说明 |
|:---:|---|---|---|---|
| **13** | `acid`, `diacid`, `polycarboxylic` | -oic acid / -dioic acid | -酸 / -二酸 | P-65.1: 羧酸为最高优先级链状 FG |
| **13** | `benzoic`, `furancarboxylic`, `pyridinecarboxylic` 等 | benzoic acid / ... | 苯甲酸 / ... | 芳环/杂环羧酸保留名，同属 carboxylic acid 类 |
| **12** | `anhydride` | -oic anhydride | -酸酐 | 羧酸酐，略低于 free acid |
| **12** | `sulfonic_acid`, `boronic` | sulfonic acid / boronic acid | 磺酸 / 硼酸 | S/B 酸类，优先级等同 anhydride |
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
| **2** | `ether`, `sulfide`, `sulfoxide` | ether / sulfide / sulfoxide | 醚 / 硫醚 / 亚砜 | 最低优先级——始终为前缀，除非无更高 FG |
| **0** | `alkane`, `alkene`, `alkyne`, `cycloalkane`, `bridged` | -ane / -ene / -yne | -烷 / -烯 / -炔 | 无 principal FG 的烃类母体 |

> **源:** `src/namepredict/layer2/kind_registry.py:21-36` (链状 FG), `37-49` (芳烃命名 FG), `50-57` (五元杂环羧酸), `71-87` (杂项环 FG)

---

## KindMeta 数据结构

`kind_registry.py:7-16` 定义了 `KindMeta` 数据类，作为所有母体种类的统一元数据容器：

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

所有母体种类通过模块级 `_bootstrap()` 函数统一注册（`kind_registry.py:314-327`），注册顺序为链状 FG → 芳烃 FG → 杂环 FG → 环烷 FG → bridged → 饱和杂环 → ScaffoldSpec。最后加载的 `_load_from_scaffold_specs()` 来自 `scaffold/specs.py`，是 **stem 的最终权威来源**——如果已有注册条目，ScaffoldSpec 会覆盖。

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

母体候选的生成当前以 **P-44 规则驱动管线**为主（`rule_driven_parent_candidates`，见
[[architecture/layer2-parent-selector]]），评分仍由 `scoring.py` 的 11 维 tuple 仲裁。经典的三类注册表
中，`fg_try_fns()` 已随 fg 注册层删除；`ring_try_fns()` / `unsat_try_fns()` 仍被保留（前者供骨架
scaffold 标注复用）。

> **源:** `src/namepredict/layer2/scoring.py:73-95`, `kind_registry.py:175-227` | 详情见 [[architecture/layer2-parent-selector]]

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

历史上共享的 BAD 元组常量（`_CORE_BAD` / `_DIACID_BAD` / `_DIONE_BAD` 等，定义于原 `parent_selector_common.py`）**已随重构删除**——各候选的互斥 keys 现在由调用方**内联**传入。典型用法：

```python
# parent_core.py —— 不饱和 FG 候选
def _ok_unsat_fg(info, flag, ekey, bad):
    return _open_chain_unsat_atoms(mol, c, db) and _no_fgs(info, bad)
```

`_no_fgs` 是唯一共享的互斥原语；每类母体的排斥集合内联在各自 producer 的 `bad` 参数里（如 `parent_core.py` / `alkynoic.py`）。这是一种声明式的约束表达：**`keys` 元组定义了"该母体不容忍的 FG 集合"**。

> **源:** `src/namepredict/layer2/fg_helpers.py`, `src/namepredict/layer2/parent_core.py`, `src/namepredict/layer2/alkynoic.py`

---

## 注册体系：kind_registry 作为单一权威来源

`src/namepredict/layer2/kind_registry.py` 是整个 FG 元数据系统的**单一权威注册中心（Single Registry Authority）**。所有对 FG 优先级的查询必须通过此模块的公共 API 进行：

| API | 功能 |
|---|---|
| `kind_registry.fg_rank(kind)` | 返回 kind 的 `fg_rank` 值 |
| `kind_registry.has_principal_fg(kind)` | 是否有 principal FG（`fg_rank > 0`） |
| `kind_registry.is_hetero_ring(kind)` / `is_carbo_ring(kind)` | 环系类型查询 |
| `kind_registry.retained_bonus(kind)` | 是否为保留名 |
| `kind_registry.parent_names(kind)` | 返回 `(en_stem, zh_stem)` 元组 |
| `kind_registry.ring_try_fns()` / `unsat_try_fns()` | 已删除——try 函数注册表整体移除 |

新增 FG 种类时，若走 principal typed 管线（acid/alcohol/ketone/amine 等），在 `principal_expression.py`
的 `_RETAINED_RING_KINDS` / `_CHAIN_KINDS` 等表达表中声明 kind 与字段即可；经典 builder
（`_open_chain_expression` / `_special_expression`）已随注册层删除。评分、候选收集与组装逻辑不变，
仍保持**开放-封闭原则**。

> **源:** `src/namepredict/layer2/kind_registry.py:89-227`

---

## 总结

官能团优先级是 NamePredict 命名正确性的基石。从 Layer 1 的排他性检测、Layer 2 的 `fg_rank` 评分、Layer 5 的后缀分派，到互斥检查规则，整个系统严格遵循 IUPAC P-41 和 P-44 的层级化思维：

1. **检测层保证"一碳一 FG"**——排他性模式使得同一碳原子不会被重复识别为多个 FG
2. **评分层保证"优先级强制"**——`fg_rank` 在 P-44 tuple 中高于环数和碳链长度，确保羧酸绝对优先于酮
3. **组装层保证"劣后 FG 降级"**——低优先级 FG 从候选后缀降为取代基前缀
4. **互斥层保证"互斥正确"**——无法共存的 FG 组合通过 `_no_fgs` 互斥检查（`fg_helpers.py`）在候选阶段即被排除
