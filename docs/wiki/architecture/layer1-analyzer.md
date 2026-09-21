# Layer1: 官能团分析 (Functional Group Analyzer)

> **源文件:** `layer1/` 目录 6 个 `.py`（793 行） | **对外接口:** `analyze`、`oxoacid_entries`、`inventory_from_info` | **下游:** L2 母体选择、L3 取代基、L4 位次、L5 词尾

## 概述

Layer1 是 SMILES → IUPAC 双语命名管线的第一层：输入 RDKit `Mol`，输出分析结果 dict（11 键）。职责三块：官能团事实（局部环境 SMARTS 命中 → 非局部后处理定型 → P-41 仲裁 → `FunctionalGroupInventory`）、结构事实（C=C / C≡C 键清单、碳索引、SSSR 环与环系拓扑）、跨层元数据（`FgSpec` 声明每个 FG 的 P-41 等级、P-43 路径、表达类型与锚点字段）。

数据流：`analyze` → `_info`（碳索引 + `_collect_fgs` + `_ring_meta`）；`_collect_fgs` 由 `_cc_bond_entries` 出键清单、`_arbitrate_parts(_detect_parts(mol), mol)` 出条目与降级集，再由 `build_inventory` 装箱。`_detect_parts` 先算 acyl 头、再铺 `_local_entries`、`radical` / `aldehyde` 去头、并入 `oxoacid_lists`、`_drop_claimed_cations` 收尾。

| 接口 | 位置 | 说明 |
| --- | --- | --- |
| `analyze` / `oxoacid_entries` / `oxoacid_lists` / `boronic_entries` | `analyzer` | 唯一入口；含氧酸条目、按 P-41 类键归位、硼酸扫描 |
| `FG_SMARTS` / `match_local_fg` / `FG_SPECS` / `oxoacid_is_acid` | `fg_local_smarts`、`fg_registry` | 局部检测表与匹配；跨层元数据与酸式判据 |
| `inventory_from_info` / `OXO_FG_CLASSES` | `functional_group_inventory` | 清单出口、含氧酸类集合 |
| `sssr_rings` / `build_ring_systems` / `kekulized` | `ring_systems` | 各层共用环访问器 |

职责边界：本层只回答"分子里有什么官能团、各 FG 归属哪些原子"，不选母体、不定位次、不生成词。SMARTS 与优先级常量均为模块级数据表，扩展官能团只加表项，检测代码不动；`__init__` 仅声明包入口。

## 核心逻辑

### SMARTS 检测表

`fg_local_smarts.FG_SMARTS` 是 `(FG 键, SMARTS)` 元组表（唯一事实来源），模式首原子即中心原子，FG 键与 `fg_registry.FG_SPECS` 一致。`_compile()` 在导入时整表编译（SMARTS 写错立即抛 `ValueError`）并归并为 `_PATTERNS_BY_FG`；`match_local_fg` 对每个 FG 取同名模式并集，按中心原子去重后升序返回 `{FG 键: [(中心原子, …), …]}`。

共享臂词常量供模式拼接：`_ACID_O`（酸式氧 = 单碳羟基氧或羧酸根阴离子氧）、`_NOT_ACID`、`_NOT_ACYCLIC_ESTER`、`_NOT_HALO`、`_ONE_C`、`_RING_HET`。

全表 35 条模式、14 个 FG 键（`sulfonamide` 由含氧酸条目派生，不出现在本表）：

| FG 键 | 条数 | 模式要点 |
| --- | --- | --- |
| radical | 1 | 哑原子 `[#0]` 的重原子邻居 |
| acyl | 1 | 哑原子所连的酰基头羰基碳 |
| acid | 1 | 羰基 + 酸式氧（羧酸） |
| ester | 2 | 羧酸酯（非环烷氧基氧）；硫代羧酸 S-酯（酯氧换 S，P-65.6.3.3.7.1） |
| acyl_halide | 1 | 羰基连卤素，排除非环酯 |
| amide | 1 | 非环 N，除羰基碳/H 外只容 C 或不连碳的 O |
| aldehyde | 1 | 至多一个碳邻居、无卤素 |
| nitrile | 1 | `+0` 的 C≡N（两原子皆限中性） |
| ketone | 3 | 双碳；单碳须连环内杂原子（P-66.1.1）；环内零碳排非内酯型酯 |
| oxoacid | 10 | 中心 P/S 含氧酸，见下节 |
| alcohol / thiol | 各 1 | 羟基连非酰基碳 / 巯基连碳 |
| amine | 3 | 按取代度分三条（H2,H3,H4 / H1 / H0），中心限 `+0`、环内与芳香氮不入 |
| cation | 8 | `_cation_local(z)` 对 `CATION_Z` 逐元素一条，见下节 |

含氧酸 10 条覆盖 P 的磷酸 P(=O)(O)₃ 与膦酸 P(=O)(O)₂-X，S 的磺酸族（磺酸/磺酸盐、磺酰卤、磺酰胺、磺酸酯）与硫酸族（硫酸氢酯/硫酸酯、硫酸/硫酸根、二硫酸链端节、多硫酸链中节）。臂词：`_O_PHOS`（单键氧三态 OH / O⁻ / O-R，臂根限碳或磷，故 P-O-P 桥氧即多聚磷酸）、`_O_ARM`、`_O_S_ARM`（S-O-S 桥臂）、`_HALO_ARM`、`_N_ARM`、`_C_ARM`、`_SNY_ARM`（膦酸第三臂）。模式只判"中心元素 + 双键氧数 + 臂型"，具体种类交由后处理归一。

阳离子 8 条由 `CATION_Z`（N / P / O / S 与四卤，P-73.1.1.1 表 7.3 的 15/16/17 族）逐元素生成；`_cation_local(z)` 的局部模式限中心 `+1`、非环，且半径 1 内无负形式电荷原子——该闸排除硝基/叠氮/N-氧化物/异氰/高氯酸根等阳离子寄居他基团的结构。

### 非局部后处理：含氧酸合一管线

检测 — 归一 — 分类三段共用一条管线，中心元素 B 与 P、S 走同一张表。

- `_oxo_partition(mol, z_idx)` 把中心原子邻居按角色分桶：双键氧 `oxo`、羟基氧 `oh`、阴离子氧 `om`、O-臂氧 `o_arm` 与臂根碳 `o_arm_root`、碳 `c`、卤素 `hal`、氮 `n`、硫 `s`（`o_arm` 存氧、`o_arm_root` 存臂根碳）。
- `_oxo_kind(mol, z_idx, part)` 依元素与角色归一 kind。B 须恰一个碳臂且两氧彻底酸式（`oh + om == 2`，硼酸酯的 O-臂不在此列）方为 `boronic`；P 由 `_OXO_KIND_P` 按 `(双键氧数, 有无直连碳臂)` 映射 `phosphate` / `phosphonate`；S 须两个双键氧，无直连碳时要求「酸式氧 + O-臂 + S 桥臂合计为 2」判 `sulfate`，有直连碳时另一臂经 `_OXO_KIND_S_ARM` 映射 `halo`→`sulfonyl_chloride`、`n`→`sulfonamide`、`acid`→`sulfonic`、`o_arm`→`sulfonate`。
- `_arm_single_attach(mol, roots, core)` 要求各臂单点回接且互不相连（`_arm_component` 做不穿 core 的连通搜索）。
- `_oxoacid_entry(mol, core)` 组装条目并归一 `oxo_kind`。中心自任母体的 `constants.OXO_CENTER_KINDS`（`phosphate` / `phosphonate` / `sulfate` / `boronic`）额外要求整分子纯度：全部重原子 = core ∪ 臂（排除臂间成环），`center_idx` 即中心原子；碳锚定的 kind 要求直连碳臂唯一并以该碳为 `center_idx`，`surr_idx` 只收阴离子氧供 L2 补 anion 标志。
- `oxoacid_entries` 按 `oxo_z` 去重排序，对缩合含氧酸（P-O-P / S-O-S）两级门控：链上须留有全酸式末端（同一中心 ≥2 个酸式氧，P-67.2.1）；无此末端时带桥氧的链节只在酸式氧（`n_oh + n_om`）达全链最大值时保留作母体，否则母体落到甲烷；链内中心让位于更少质子化的酸中心（`_oxo_bridge_arms` 计单键 O 上另一同元素中心数，非 P/S 中心恒 0）。保留的链节即送入 [[architecture/layer2-parent-selector]] 的母体候选。
- `boronic_entries` 按元素直接扫 B（B 无局部双键氧 SMARTS 可匹配），以 `(B, 其氧邻居)` 为 core 走同一 `_oxoacid_entry`；`oxoacid_lists` 合并两者，按 `_OXO_CLASS_BY_KIND` 归位到 P-41 类别键：磺酰胺落 `"sulfonamide"`（类 11），其余落 `"oxoacid"`（类 9）。
- `_PRESENCE_SKIP` 含 `oxoacid` 与 `sulfonamide`：含氧酸不参与存在性判定，纳入会改写压制结果。

### 非局部后处理：其余 FG

- `_detect_parts` 里 acyl 头先算：头碳集 `heads` 从 `radical` 与 `aldehyde` 候选中剔除（哑原子标记的酰基位点不重复计）。
- `_local_entries` 按 `_LOCAL_ENTRY_FGS` 顺序（acid / alcohol / ester / amide / ketone / amine / thiol / nitrile / acyl_halide / cation）为纯局部 FG 组装条目；`_fg_entry` / `_surr_idx` 给出中心与重原子周边，碳中心不在环内时排除环内邻居；`_radical_entry` 无周边。`_drop_claimed_cations` 剔除中心已被其它 FG 检测器命中的阳离子条目。
- `_cc_bond_entries` 一遍遍历全部键，按 C=C（排除芳香键）与 C≡C 分别组装条目，`{"c1","c2"}` 为有序碳对。
- `_arbitrate_parts(parts, mol)` 做 P-41 仲裁：`present` 为有条目且不在 `_PRESENCE_SKIP` 的 FG；分子含负形式电荷原子（`_has_negative_atom`）时从 `present` 剔除 `cation`（阴离子类 4 > 阳离子类 6）。`_SUPPRESSIBLE`（acid / ester / acyl_halide / amide / nitrile / aldehyde）遇更高优先级 FG 即退出：`_LEAF_DEMOTED`（acid / nitrile）整组碳排除出主链并标 `demoted`（P-61.1.3 carboxy / cyano），其余（酯/酰胺/醛/酰卤）置空列表降为"氧代"类前缀，由 L3 锚定叶识别。

### 清单与 occurrence

`build_inventory(lists, mol, demoted)` 按 `_FG_KEYS`（`FG_SPECS` 的 `fg` 投影，产出顺序即注册顺序）逐条装箱为 `FunctionalGroupInventory`：`entries` 为 `FunctionalGroupOccurrence`（`id` = `"类别:序号"`、`group_class`、`characteristic_atoms`、`parent_anchors`、`payload`、`demoted`），另带 `has_anion`。`occurrences(group_class)` 按类查询且不含降级条目，`demoted_entries()` 取全部降级叶。

特征原子默认取 `center_surr_atoms`（中心 ∪ `surr_idx`）；`FG_ATOM_FNS` 登记例外：`oxoacid` 与 `sulfonamide` 共用 `_oxoacid_atoms`（中心 + 锚点碳 + 中心的非碳邻居，碳臂留给链/取代基侧），`cation` 用 `_cation_atoms`（只占中心一个原子）。`parent_anchors` 由 `_ANCHOR_KEYS`（`FG_SPECS.anchors` 投影）经 `_indices` 收集。`OXO_FG_CLASSES` = { `OXOACID`, `SULFONAMIDE` }，即含氧酸合一类的全部 P-41 类别。

`inventory_from_info(info)` 是下游唯一入口，缺 `fg_inventory` 时抛 `KeyError`。参见 [[concepts/atom-ownership]]。

### 环系拓扑

`_ring_meta` 汇总 `rings`、`n_rings`、`has_ring`、`ring_systems`、`n_ring_systems`。`build_ring_systems(mol)` 以 `sssr_rings` 为输入，`_ring_pairs` 单遍扫环对（共享 ≥2 原子入 `fusion_edges`、恰好 1 个入螺环对），`_components` 以 `fusion_edges + spiro_edges` 建邻接表求连通分量（分量按最小成员升序），螺环对因此并入同一环系。`_system_entry` 为一分量组装环系 dict：`atom_ids`（环成员原子，升序）、`sssr_indices`（SSSR 环索引）、`fusion_edges`（`(i, j, 共享原子升序)`）、`spiro_edges`（`(i, j, 螺原子)`）、`free_spiro_atoms`（自由螺原子，P-24.1，升序）、`n_rings` / `n_atoms`（环数 / 环系原子数）、`hetero_atoms`（环系内非碳非氢原子 `{"idx","Z"}`，升序）、`topology`（有自由螺原子时为 `"spiro"`，否则 `None`）。

自由螺原子由 `_free_spiro_atoms` 判定：取分量环子图的邻接表，`_split_count` 数去掉某螺原子后其环内邻居分属的连通分量数（只沿环内邻居走、不穿 `removed`），≥2 即该原子去掉后环子图断成两块。

`sssr_rings` 是各层统一环访问器，按分子经 `tools.memo.by_mol` 记忆；`kekulized` 返回去芳香标志的 Kekulé 副本，失败返 `None`。参见 [[guides/adding-new-ring-system]]。

## 与其他层的契约

- `analyze` 的 11 键（`mol`、`carbon_ids`、`n_carbons`、`double_bonds`、`triple_bonds`、`fg_inventory`、`rings`、`n_rings`、`has_ring`、`ring_systems`、`n_ring_systems`）是 L2–L5 的唯一事实来源；官能团事实只经 `fg_inventory` 出口，`has_ring` 被 `namer` 用作链/环分派。全表见 [[reference/core-data-contracts]]。
- `double_bonds` / `triple_bonds` 的 `{"c1","c2"}` 有序碳对由 L2 `principal_expression._unsat_bond_fields` 消费为 `double_bond` / `double_bonds` 等字段。
- 环系 dict 的 `sssr_indices` / `fusion_edges` / `spiro_edges` / `free_spiro_atoms` / `atom_ids` 是 L2 拆解（`fused_system` / `spiro_system` / `bridged_system`）、L4 编号与 L5 环名词的输入。
- `FgSpec` 7 字段（`fg` / `p41` / `path` / `expr` / `anchors` / `parent_anchor_fields` / `locant_source`）被 L2 优先级、L4 位次（`locant_calc` 由 `FG_SPECS` 投影）、L5 前缀派生；新增类别见 [[guides/adding-new-functional-group]]。
- 酸式判据 `oxoacid_is_acid`（配 `OXO_ACID_P41` / `OXO_ACID_KINDS` / `OXO_ACID_KIND_BY_H`）由 L2 `principal._effective_priority` 消费：候选为 `oxoacid` 且全部 occurrence 均酸式时升到类 7（优先于类 9 酯与类 11 酰胺）。
- `constants.OXO_CENTER_KINDS` 判"中心自任母体"：L2 `parent_select` 只保留 FG 原子、L5 `assembler` 走功能母体词尾，`boronic` 的 `n_oh` 决定硼酸/硼酸氢/硼酸酯词尾。
- `cation` 类（类 6）与 `has_anion` 供 L2 `principal` 判阴离子让位；L2 `parent_skeleton` 对全为阳离子的主基团只向阳离子原子收敛。
- `OXO_FG_CLASSES` 供 L2 `principal_expression` 判阴离子标志、`parent_select` 归并特征原子。
- `match_local_fg`、`oxoacid_entries` 与私有工具 `_oxo_bridge_arms`、`_double_bonded_o_idxs`、`_alkoxy_c_of` 被 L2 `parent_select`、`principal_expression` 跨层复用。
