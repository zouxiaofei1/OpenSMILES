# 原子归属 (Atom Ownership)

> **概念层级:** 跨层核心概念 | **涉及层次:** Layer2 (设置边界)、Layer3 (消耗边界)
> **核心数据结构:** `owned_atoms: frozenset[int]` | `ClaimedBlock` | `CoverageLedger`

---

## 概述

在 NamePredict 的 6 层命名流水线中，**原子归属 (atom ownership)** 是一个贯穿 Layer2 和 Layer3 的基础概念。它定义了分子中每个重原子（即非氢原子）的"所有权"——哪些原子属于母体 (parent)，哪些原子属于取代基 (substituent)，以及是否存在 gap（未被任何实体声明的原子）或 overlap（被多个实体重复声明的原子）。

IUPAC 命名的本质是：选定一个母体结构作为骨架，然后将分子其余部分描述为连接在该骨架上的取代基。原子归属机制精确地划定了这条边界，确保每个原子被恰好一次 ——不遗漏、不重复——分配给了正确的命名实体。

原子归属由 **Layer2 设置边界**，**Layer3 在边界外提取取代基**，**Coverage Ledger 在 Layer3 输出后进行最终验证**。

> **源:** `src/namepredict/layer2/parent_ownership.py:279-282`
> ```python
> def compute_owned_atoms(parent: dict, mol: Mol) -> frozenset[int]:
>     """Union chain + kind-specific FG atoms (terminal ownership set)."""
>     return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))
> ```

## 母体拥有的原子 (Owned Atoms)

### 计算模型

母体拥有的原子集合 `owned_atoms` 是两个子集的 **并集**：

1. **骨架链原子 (Chain Atoms):** 母体主链上的所有碳原子索引，由 `parent["chain"]` 字段提供。这是母体的碳骨架。
2. **官能团异原子 (FG Heteroatoms):** 母体主要官能团中非碳的重原子（O、N、S、Cl 等），由 `_kind_fg_atoms()` 函数根据母体类型 `parent["kind"]` 动态计算。

计算公式：
```
owned_atoms = chain_atoms ∪ fg_atoms
```

> **源:** `src/namepredict/layer2/parent_ownership.py:279-282`

### 不可变性

`owned_atoms` 被存储为 Python `frozenset[int]` 类型，在 Layer2 最终确定后 **不可修改**。这个不可变性保证了 Layer3 提取取代基时不会意外地修改母体边界，也避免了官能团原子的重复计入。

> **源:** `src/namepredict/layer2/parent_ownership.py:284-287`
> ```python
> def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
>     if isinstance(parent.get("owned_atoms"), frozenset):
>         return parent
>     return {**parent, "owned_atoms": compute_owned_atoms(parent, mol)}
> ```

若 `parent` 字典中已经存在 `frozenset` 类型的 `owned_atoms` 字段，则直接返回原字典（幂等保障），避免重复计算。

### 骨架链原子

作为所有母体的共同基础，骨架链原子由 `_chain_atoms()` 函数获取：

> **源:** `src/namepredict/layer2/parent_ownership.py:9`
> ```python
> def _chain_atoms(parent: dict) -> set[int]:
>     return set(parent.get("chain") or [])
> ```

对于链状母体（烷烃、烯烃、炔烃等），`chain` 是主链碳原子列表；对于环状母体，`chain` 可能为空（环原子通过 `ring` 字段管理），此时骨架原子通过环相关的逻辑单独纳入。

### 按 FG 类型的异原子归属

不同类型的母体主要官能团拥有不同的异原子集合。以下按母体类型说明各自的归属范围：

#### 羧酸 (Carboxylic Acid)

羧酸母体拥有 **每个羧基碳 + C=O 双键氧 + O-H 单键氧**。支持单羧酸（`cooh_c_idx`）和多羧酸（`cooh_c_idxs`）两种情况。

> **源:** `src/namepredict/layer2/parent_ownership.py:45-68`

实际归属: `{C_carboxyl} ∪ {=O} ∪ {-OH}` 对于每个羧基。`_acid_o_atoms()` 通过检测碳原子的双键氧邻居 (`_dbl_o_idx`) 和单键氧邻居 (`_single_o_idx`) 来确定这两个氧原子，`_acid_fg_atoms()` 遍历 `_cooh_c_idxs()` 合并所有羧基。

#### 醛 (Aldehyde)

醛母体拥有 **醛基碳 + C=O 双键氧**。支持单醛（`aldehyde_c_idx`）与多醛（`aldehyde_c_idxs`）——环上外环 -CHO 的多醛（-dicarbaldehyde）也取全部醛基碳。

> **源:** `src/namepredict/layer2/parent_ownership.py:70-82`

归属（每个醛基）: `{C_aldehyde} ∪ {=O}`。

#### 酮 (Ketone)

酮母体拥有 **羰基碳 + C=O 双键氧**。对于乙酰基类酮（如 acetophenone），还会额外包含 **乙酰甲基碳** (`acetyl_methyl_idx`)。

> **源:** `src/namepredict/layer2/parent_ownership.py:130-140`

归属: `{C_carbonyl} ∪ {=O} ∪ {acetyl_methyl}`（如存在）。

#### 醇/酚 (Alcohol/Phenol)

醇类母体拥有 **每个连接碳 + 对应的羟基氧**。支持单醇（`oh_c_idx`）和多元醇（`oh_c_idxs`）。

> **源:** `src/namepredict/layer2/parent_ownership.py:112-128`

归属（每个 OH 基团）: `{C_attach} ∪ {-OH}`。

#### 胺/苯胺 (Amine/Aniline)

胺类母体拥有 **每个连接碳 + 与其单个键合的氮原子**。支持单胺和多胺。

> **源:** `src/namepredict/layer2/parent_ownership.py:142-157`

归属: `{C_attach} ∪ {N}`。注意：氮上连接的氢原子不算重原子，不计入 `owned_atoms`；氮上若连接额外的碳链取代基，那些碳在母体边界之外，将由 Layer3 提取。

#### 酰胺 (Amide)

酰胺母体拥有 **酰胺羰基碳 + C=O 双键氧 + 酰胺氮**。

> **源:** `src/namepredict/layer2/parent_ownership.py:30-38`

归属: `{C_amide} ∪ {=O} ∪ {N}`。`_amide_n_from_c()` 通过查找羰基碳上原子序数为 7 的邻居来确定酰胺氮。

#### 酯 (Ester)

酯母体拥有 **酯羰基碳 + C=O 双键氧 + 酯氧 (-O-)**。当前实现只覆盖单酯/苯甲酸酯（`ester_c_idx`），不纳入烷氧臂原子。

> **源:** `src/namepredict/layer2/parent_ownership.py:185-208`

归属: `{C_ester} ∪ {=O} ∪ {-O-}`。`_one_ester_fg` 加入羰基碳、双键氧与单键酯氧（acid 侧，不含 alkoxy 臂）。

#### 醚 (Ether)

醚母体拥有 **醚氧 + 两侧碳臂完整碳链**。两侧碳臂均通过 `_longest_from` 走最长路径获取。

> **源:** `src/namepredict/layer2/parent_ownership.py:84-105`

归属: `{O_ether} ∪ {arm1_atoms} ∪ {arm2_atoms}`。

#### 硫醚 (Sulfide)

硫醚母体的归属模式与醚完全对称，将 O 替换为 S：`{S_sulfide} ∪ {arm1_atoms} ∪ {arm2_atoms}`。

> **源:** `src/namepredict/layer2/parent_ownership.py:90-110`

#### 硫醇 (Thiol)

硫醇母体拥有 **连接碳 + 巯基硫**。

> **源:** `src/namepredict/layer2/parent_ownership.py:171-183`

归属: `{C_attach} ∪ {S}`。

#### 腈 (Nitrile)

腈母体拥有 **氰基碳 + 三键氮**。

> **源:** `src/namepredict/layer2/parent_ownership.py:159-169`

归属: `{C_nitrile} ∪ {N≡}`。

#### 酸酐 (Anhydride)

酸酐母体拥有 **两个酰基碳 + 桥氧 + 两个羰基氧**。

> **源:** `src/namepredict/layer2/parent_ownership.py:210-220`

归属: `{C_acyl1} ∪ {C_acyl2} ∪ {=O1} ∪ {=O2} ∪ {O_bridge}`。

#### 酰卤 (Acyl Halide)

酰卤母体拥有 **酰基碳 + 羰基氧 + 卤素 (F/Cl/Br/I)**。通过 `cl_idx` 或 `hal_idx` 字段定位卤原子。

> **源:** `src/namepredict/layer2/parent_ownership.py`

归属: `{C_acyl} ∪ {=O} ∪ {X}`。

#### 酰基残基 (Acyl)

`*`-锚定酰基残基/环外酰基头（苯甲酰、furan-2-carbonyl，`_acyl_fg_atoms` `parent_ownership.py:232`）母体拥有 **羰基头碳 + 羰基氧**（=O 归母体，不落入 oxo 前缀）。

> **源:** `src/namepredict/layer2/parent_ownership.py:232`

归属: `{C_acyl_head} ∪ {=O}`。

#### 磷酸 (Phosphate)

磷酸/磷酸酯母体（`kind == "phosphate"`）拥有 **P 中心 + 全部 4 个氧**（1 个 =O 与 3 个单键 O，含 O–R 桥氧）。桥氧归母体，O 上的烷基臂在母体边界之外，由 Layer3 以 `o_side` 取代基提取（`claim_extract._ESTER_O_SIDE_KINDS` 含 `"phosphate"`），最终由 L5 `phosphate_names` 拼进磷酸酯整名。

> **源:** `src/namepredict/layer2/parent_ownership.py:242-252`

归属: `{P} ∪ {=O} ∪ {O_single × 3}`。

### FG 原子聚合

所有 FG 特定的原子集合通过 `_kind_fg_atoms()` 函数统一聚合并集。它的实现是纯"字段驱动"的——检查 parent 字典中是否存在某个 FG 定位字段，若存在则调用对应的子函数计算原子集，最后取所有非空子集的并集：

> **源:** `src/namepredict/layer2/parent_ownership.py:254-277`

> 注：12 个扩展 FG（磺酸/亚砜/砜/硼酸/氨基甲酸酯/脲/胍/肼等）未实现，其归属函数（`_sulfonic_*`/`_sulfoxide_*`/`_boronic_*`/`_carbamate_*`/`_urea_*` 等）不存在；磷酸已实现（`_phosphate_fg_atoms`，`parent_ownership.py:242`）。

这种字段驱动设计使得添加新的 FG 类型只需新增两个函数（FG 原子计算函数 + 条件判断），不影响其他类型的逻辑。

## ClaimedBlock 与 SideSlot

### 从母体边界到取代基声明

当 `owned_atoms` 确定后，Layer3 的工作是将所有 **不在** `owned_atoms` 中的重原子组件识别为取代基。`ClaimedBlock` 是这一过程的中间数据结构，它记录了取代基原子块的基本拓扑信息。

### ClaimedBlock 数据结构

```python
@dataclass(frozen=True)
class ClaimedBlock:
    slot: SideSlot       # 母体侧连接点的化学角色
    attach_parent: int   # 母体上被连接的原子索引
    root: int            # 取代基侧的根原子索引
    atoms: frozenset[int]  # 取代基包含的所有重原子
```

> **源:** `src/namepredict/layer3/claimable_block.py:23-28`

`ClaimedBlock` 是一个冻结的 dataclass，不可变，保证在命名过程中不会被意外修改。

### SideSlot 枚举

`SideSlot` 描述了母体侧连接点的化学环境，影响取代基的命名优先级和排序：

| 枚举值 | 含义 | 判断条件 |
|--------|------|----------|
| `CHAIN_C` | 母体链碳上的连接 | 连接原子是碳、不在环中、非酰胺/胺/醚角色 |
| `RING_C` | 母体环碳上的连接 | 连接原子是碳、在环中 |
| `AMIDE_N` | 酰胺氮上的连接 | 连接原子是 N，单键连接到母体中含有双键氧的羰基碳 |
| `AMINE_N` | 胺氮上的连接（8584795 新增） | 连接原子是非芳香 N，至少一个邻居是母体内的非羰基碳 |
| `ETHER_O` | 醚氧上的连接 (已被取代基臂占据) | 连接原子是 O，恰好两个碳邻居（无氢邻居） |
| `OTHER` | 其他杂原子的连接 | 不满足以上任何条件 |

> **源:** `src/namepredict/layer3/claimable_block.py:12-20, 84-96`

SideSlot 的推导由 `derive_slot()` 函数完成，它仅根据母体侧已归属原子的化学环境判断槽位类型。例如，检测酰胺氮的逻辑 (`_is_amide_n`) 会验证：该原子为氮 (Z=7)、单键连接到一个已归属的羰基碳，且该碳带有双键氧。`AMINE_N` 的推导 (`_is_amine_n`) 验证：非芳香 N 至少一个邻居是 owned 内非羰基碳——供 L5 生成 `N-` 前缀取代基（P-62.2）。

### 声明块的生成

`claim_block()` 函数负责创建单个 ClaimedBlock：

1. 验证 `attach_parent` 确实在 `owned_atoms` 内（否则连接点不在母体中，无效）
2. 验证 `root` 确实在母体外且是重原子 (`_is_outside_root`)
3. 通过 `cut_block()` 将整个外部组件从 root 开始切割出来
4. 验证该组件仅有 **唯一一个** 连接回母体的点（多连接点=桥连=不是简单取代基，返回 None）

> **源:** `src/namepredict/layer3/claimable_block.py:116-133`

然后 `iter_claims()` 遍历所有外部重原子组件（通过 `side_roots()` 获取所有母体外的起点，`cut_block()` 切割出组件，`_unique_components()` 去重），并为每个组件创建 ClaimedBlock。结果按 `(attach_parent, root, slot)` 三元组排序，保证输出顺序确定。

> **源:** `src/namepredict/layer3/claimable_block.py:189-196`

## Coverage Ledger: Gap 与 Overlap 验证

### 完整性约束

原子归属的正确性由一个核心约束定义：**分子的每个重原子必须恰好被一个命名实体覆盖**。这意味着：

- 母体声明了其骨架链 + 官能团异原子的所有权
- 每个取代基声明了其原子块的所有权
- 两个集合的并集必须覆盖所有重原子（无 gap）
- 两个集合的交集必须为空（无 overlap）

### CoverageLedger 数据结构

```python
@dataclass(frozen=True)
class CoverageLedger:
    owned_atoms: frozenset[int]
    named_claims: tuple[SubstituentName, ...]
    gap: frozenset[int]       # 未被任何实体覆盖的重原子
    overlap: frozenset[int]   # 被多个实体重复覆盖的重原子

    @property
    def complete(self) -> bool:
        return not self.gap and not self.overlap
```

> **源:** `src/namepredict/layer3/coverage.py:13-23`

### 构建过程

`build_coverage_ledger()` 函数构造覆盖台账：

1. 从 RDKit Mol 对象获取所有重原子集合 (`_heavy_atoms`)
2. 计算所有被覆盖原子的并集：`owned_atoms ∪ (每个 substituent 的 claim.atoms)`
3. gap = 所有重原子 - 被覆盖原子
4. 计算每个原子的"被计数次数"（在 owned_atoms 中计 1 次，在每个取代基中计对应的次数）
5. overlap = 计数次数 > 1 的原子集合

> **源:** `src/namepredict/layer3/coverage.py:57-72`

### gap 与 overlap 的含义

- **gap (缺口):** 存在某个重原子既不属于母体也不属于任何取代基。这通常意味着某个官能团或侧链未能被正确提取，或者母体边界定义不完整。gap 原子是 Layer3 提取逻辑的 bug 信号。

- **overlap (重叠):** 存在某个重原子同时被母体和某个取代基声明，或被多个取代基重复声明。这通常意味着：
  - Layer2 的 FG 异原子归属有遗漏，导致某个原子未能被母体声明，反而被取代基错误地提取
  - 或者取代基切割逻辑未能正确地尊重母体边界
  - 或者存在功能性基团冲突（如母体是醇，但链上另一个 OH 被同时视为母体的一部分和取代基）

### 在流水线中的位置

Coverage Ledger 在 Layer3 提取取代基之后构建，在 Layer4（编号）之前验证。在 `namer.py` 中：

```python
subst = extract_substituents(info, parent, name_mode=name_mode)
complete = _ledger_complete(mol, parent["owned_atoms"], subst)
```

> **源:** `src/namepredict/namer.py:152-153`

只有当 `CoverageLedger.complete` 为 `True` 时（无 gap、无 overlap），命名流程才会继续进入编号和组装阶段。

## 原子归属的生命周期

总结原子归属在整个流水线中的生命周期：

```
Layer1 (analyze)
  └→ 识别所有 FG，但尚未确立归属关系

Layer2 (parent_selector)
  └→ 选定母体结构
  └→ finalize_parent_ownership() 计算 owned_atoms
      ├── _chain_atoms: 骨架碳链
      └── _kind_fg_atoms: FG 异原子（按母体类型）
  └→ owned_atoms 以 frozenset 冻结，不可变更

Layer2/3 边界 (claimable_block)
  └→ iter_claims() 遍历 owned_atoms 外的所有重原子组件
  └→ 为每个组件创建 ClaimedBlock (slot + attach_parent + root + atoms)

Layer3 (substituent_extractor)
  └→ 以 owned_atoms 为边界，在边界外提取取代基
  └→ 每个取代基的命名基于其 ClaimedBlock
  └→ build_coverage_ledger() 验证 gap/overlap
  └→ CoverageLedger.complete == True → 继续到 layer4

Layer4 → Layer5
  └→ 编号和名称组装不再涉及原子归属（边界已在之前确立）
```

## 相关页面

- [[architecture/layer2-parent-selector]] — Layer2 母体选择器，决定 owned_atoms 的内容
- [[architecture/layer3-substituents]] — Layer3 取代基提取器，基于 owned_atoms 边界提取取代基
- [[concepts/functional-group-priority]] — 官能团优先级，决定哪个 FG 成为母体主要官能团
- [[architecture/overview]] — 6 层架构总览
- [[reference/core-data-contracts]] — 核心数据契约（parent dict 字段规范）
