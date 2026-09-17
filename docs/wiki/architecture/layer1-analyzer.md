# Layer1: 官能团分析 (Functional Group Analyzer)

> **源文件:** `layer1/` 目录 6 个 `.py`（725 行） | **对外接口:** `analyze`、`oxoacid_entries`、`inventory_from_info` | **下游:** L2 母体选择、L3 取代基、L5 词尾

## 概述

Layer1 是 SMILES → IUPAC 双语命名管线的第一层：输入 RDKit `Mol`，输出一个分析结果 dict。职责有三块。

- 官能团事实：局部环境 SMARTS 命中 → 非局部后处理定型 → P-41 仲裁 → `FunctionalGroupInventory`。
- 结构事实：C=C / C≡C 键清单、碳原子索引、SSSR 环与环系拓扑（`ring_systems`）。
- 跨层元数据：`FgSpec` 声明每个 FG 的 P-41 等级、P-43 路径与锚点字段，供 L2/L4/L5 派生。

职责边界：本层只回答"分子里有什么官能团、各 FG 归属哪些原子"，不选母体、不定位次、不生成词。SMARTS 与优先级常量均为模块级数据表，扩展官能团只加表项，检测代码不动。

模块分工：`fg_local_smarts` 局部 SMARTS 表与匹配；`fg_registry` 跨层元数据与酸式判据；`analyzer` 非局部后处理、仲裁与 `analyze` 出口；`functional_group_inventory` 带类型清单与查询；`ring_systems` 环与环系拓扑（L1/L4/L5 共用）；`__init__` 仅声明包入口。

## 检测：FG_SMARTS 表驱动

`fg_local_smarts.FG_SMARTS` 是 `(FG 键, SMARTS)` 元组表（唯一事实来源），FG 键与 `fg_registry.FG_SPECS` 一致。`_compile()` 在导入时整表编译（SMARTS 写错立即抛 `ValueError`）并归并为 `_PATTERNS_BY_FG`；`match_local_fg` 对每个 FG 取同名模式并集，按中心原子（模式首原子）去重后升序返回 `{FG 键: [(中心原子, …), …]}`。

共享臂词常量供模式拼接：`_ACID_O`（酸式氧 = 单碳羟基氧或羧酸根阴离子氧）、`_NOT_ACID`、`_NOT_ACYCLIC_ESTER`、`_NOT_HALO`、`_ONE_C`、`_RING_HET`。

全表 24 条模式、13 个 FG 键（`sulfonamide` 由含氧酸条目派生，不出现在本表）：

| FG 键 | 条数 | 模式要点 |
| --- | --- | --- |
| radical | 1 | 哑原子 `[#0]` 的重原子邻居 |
| acyl | 1 | 哑原子所连的酰基头羰基碳 |
| acid | 1 | 羰基 + 酸式氧（羧酸） |
| ester | 2 | 羧酸酯（非环烷氧基氧）；硫代羧酸 S-酯（酯氧换 S，P-65.6.3.3.7.1） |
| acyl_halide | 1 | 羰基连卤素，排除非环酯 |
| amide | 1 | 非环 N，除羰基碳/H 外只容 C 或不连碳的 O |
| aldehyde | 1 | 至多一个碳邻居、无卤素 |
| nitrile | 1 | C≡N |
| ketone | 3 | 双碳；单碳须连环内杂原子（P-66.1.1）；环内零碳排非内酯型酯 |
| oxoacid | 7 | 中心 P/S 含氧酸，见下节 |
| alcohol / thiol | 各 1 | 羟基连非酰基碳 / 巯基连碳 |
| amine | 3 | 按取代度分三条（H2,H3,H4 / H1 / H0），环内与芳香氮不入 |

含氧酸 7 条覆盖：磷酸 P(=O)(O)₃、膦酸 P(=O)(O)₂-X、磺酸/磺酸盐、磺酰卤、磺酰胺、磺酸酯、硫酸氢酯/硫酸酯。臂词为 `_O_PHOS`（单键氧三态 OH / O⁻ / O-R，臂根限碳或磷，故 P-O-P 桥氧即多聚磷酸；排除哑原子臂）、`_O_ARM`、`_HALO_ARM`、`_N_ARM`（非环氮）、`_C_ARM`（直连碳）、`_SNY_ARM`（膦酸第三臂，C/S/N 均可）。模式只判"中心元素 + 双键氧数 + 臂型"，具体种类交由后处理归一。

## 非局部后处理

SMARTS 只约束中心原子邻域，整分子环境由 `analyzer` 后处理定型。

### 含氧酸：P/S 合一管线

检测 — 归一 — 分类三段共用一条管线，中心元素 P 与 S 走同一张表。

- `_oxo_roles(mol, z_idx)` 计中心原子邻居角色：双键氧 `oxo`、羟基氧 `oh`、阴离子氧 `om`、O-臂 `o_arm`、碳 `c`、卤素 `hal`、氮 `n`、硫 `s`。
- `_oxo_kind(mol, z_idx)` 依元素与角色归一 kind。P 由 `_OXO_KIND_P` 按 `(双键氧数, 有无直连碳臂)` 映射为 `phosphate` / `phosphonate`；S 须两个双键氧，无直连碳时由「一酸式一 O-R」判 `sulfate`（否则失败），有直连碳时另一臂经 `_OXO_KIND_S_ARM` 映射 `halo`→`sulfonyl_chloride`、`n`→`sulfonamide`、`acid`→`sulfonic`、`o_arm`→`sulfonate`。
- `_oxo_arm_roots` 取中心原子的 (直连碳臂, O-臂根碳, 羟基氧, 阴离子氧, 酸式氧)；`_arm_single_attach` 要求各臂单点回接且互不相连（`_arm_component` 做不穿 core 的连通搜索）。
- `_oxoacid_entry` 组装条目并归一 `oxo_kind`。中心自任母体的 `_OXO_Z_ANCHORED`（`phosphate` / `phosphonate` / `sulfate`）额外要求整分子纯度：全部重原子 = core ∪ 臂，臂间成环与焦磷酸一类均被排除，`center_idx` 即中心原子；碳锚定的 kind 要求直连碳臂唯一并以该碳为 `center_idx`，`surr_idx` 只收阴离子氧，供 L2 补 anion 标志。
- `oxoacid_entries` 按 `oxo_z` 去重排序，对缩合磷酸（P-O-P）做两级门控：链上须留有全酸式末端（同一 P ≥2 个酸式氧，P-67.2.1）；链内 P 让位于更少质子化的磷酸中心（`_p_bridge_arms` 计 O-P 桥氧数，非 P 中心恒 0）。
- `oxoacid_lists` 按 `_OXO_CLASS_BY_KIND` 把条目归位到 P-41 类别键：磺酰胺落 `"sulfonamide"`（类 11），其余落 `"oxoacid"`（类 9）。`phosphate_entries` 是「仅 `phosphate` kind」的视图。
- `_PRESENCE_SKIP` 含 `oxoacid` 与 `sulfonamide`：含氧酸不参与存在性判定，纳入会改写压制结果。

模块头注释声明：含氧酸中心的非局部判据与 `oxo_kind` 归一本模块；羰基原语 `_double_bonded_o_idxs` / `_alkoxy_c_of` 供 L2 主基团表达式复用。

### 其余 FG

- `_detect_parts` 里 acyl 头先算：头碳集 `heads` 从 `radical` 与 `aldehyde` 候选中剔除（哑原子标记的酰基位点不重复计）。
- `_local_entries` 按 `_LOCAL_ENTRY_FGS` 顺序为纯局部 FG 组装条目；`_fg_entry` / `_surr_idx` 给出中心与重原子周边，碳中心不在环内时排除环内邻居；`_radical_entry` 无周边。
- 双键与三键分别由 `_is_cc_double` / `_is_cc_triple` 过滤（排除芳香键），`_bond_entry` 输出有序碳对。

## FgSpec 跨层元数据

`fg_registry.FgSpec` 是 FG 元数据的唯一事实来源，7 个字段：`fg`（`FunctionalGroupClass` 值，同作 L1 列表 key）、`p41`（P-41 主官能团等级，0 = 非主官能团）、`path`（P-43 优先级路径）、`expr`（`suffix` / `prefix_only` / `legacy_compat`）、`anchors`（occurrence payload 锚点 key，空 = 不收集）、`parent_anchor_fields`（parent 锚点字段，单/复）、`locant_source`（`attachment` / `attachment_exocyclic` / `anchor_field`）。

`FG_SPECS` 登记 14 项，下游表由其派生：`functional_group_inventory._FG_KEYS` 取 `fg` 投影（产出顺序即注册顺序），`_ANCHOR_KEYS` 取 `anchors` 投影。要点：

- `oxoacid` 为 `p41=9 path=(0,) anchors=("center_idx",)`，中心 P/S 合一类，`oxo_kind` 由 L1 payload 归一。
- `sulfonamide` 为 `p41=11 path=(1,)`，与酰胺同组并排在其后。
- 类 7 的"酸式"判据：`OXO_ACID_P41`（=8，类 7 内排在羧酸 `acid=7` 之后）、`OXO_ACID_KINDS`（`{"sulfonic"}`）、`OXO_ACID_KIND_BY_H`（`{"phosphonate": 1}`）、判据函数 `oxoacid_is_acid(payload)`。中心带 O⁻ 按酸根处理；kind 在 `OXO_ACID_KINDS` 内按酸式；碳锚定 P 酸按 `n_oh ≥ 阈值`。
- 含氧酸为酸式时按类 7 参与主基团竞争，优先于类 9 酯与类 11 酰胺；中性磷酸/硫酸酯仍属类 9，按相应酸排羧酸酯之后。判据消费方在 L2 `principal.py`：候选类为 `oxoacid` 且其全部 occurrence 均为酸式时，实际优先级升为 `OXO_ACID_P41`。

## 统一出口

`analyze(mol)` 返回完整结果 dict：`mol`、`carbon_ids`、`n_carbons`、`double_bonds`、`triple_bonds`、`fg_inventory`、`rings`、`n_rings`、`has_ring`、`ring_systems`、`n_ring_systems`。

调用链：`analyze` → `_info`（碳索引 + `_collect_fgs` + `_ring_meta`）；`_collect_fgs` = `_arbitrate_parts(_detect_parts(mol))`，再经 `_filter_bond_entries` 收集双/三键，终由 `build_inventory` 装箱。`_detect_parts` 内先算 acyl 头、再铺 `_local_entries`、`radical` / `aldehyde` 去头，最后并入 `oxoacid_lists`。

- **P-41 仲裁**：`_arbitrate_parts` 取 `fg_registry` 的 `p41` 表，若存在更高优先级（数值更小）的 FG 则压制 `_SUPPRESSIBLE` 组合羰基与腈；`_LEAF_DEMOTED`（`acid` / `nitrile`）整组碳排除出主链并标 `demoted`，酯/酰胺/醛/酰卤则降为"氧代"类前缀，由 L3 锚定叶识别。
- **inventory**：`build_inventory(lists, mol, demoted)` 产出 `FunctionalGroupInventory`，条目为 `FunctionalGroupOccurrence`（id、`group_class`、`characteristic_atoms`、`parent_anchors`、`payload`、`demoted`）。特征原子默认取中心与周边并集；`FG_ATOM_FNS` 登记例外，`oxoacid` 与 `sulfonamide` 共用 `_oxoacid_atoms`（锚点碳 + 中心 P/S + 中心的非碳邻居，碳臂留给链/取代基侧）。`occurrences()` / `demoted_entries()` 提供按类查询。
- **取值**：`inventory_from_info(info)` 是下游唯一入口，缺 `fg_inventory` 时抛 `KeyError`。
- **环事实**：`_ring_meta` 汇总 `rings`、`n_rings`、`has_ring`、`ring_systems`、`n_ring_systems`。`ring_systems.build_ring_systems` 以 SSSR 为输入，`_ring_pairs` 单遍扫环对（共享 ≥2 原子为稠合边、恰好 1 个为螺环对），`_components` 用并查集做环索引连通分量，`_system_entry` 组装 `atom_ids`、`sssr_indices`、`fusion_edges`、`n_rings`、`n_atoms`、`hetero_atoms`（`topology` 与 `is_aromatic_mancude` 置 `None` 留待上层填充）。`sssr_rings` 是各层统一环访问器，与 `_sssr` 同一次 `tools.memo.by_mol` 记忆；`kekulized` 返回去芳香标志的 Kekulé 副本，失败返 `None`。
