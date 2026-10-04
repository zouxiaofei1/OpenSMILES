# Layer2: 母体选择 (Parent Selection)

> **源文件:** `src/opensmiles/layer2/` 12 个 `.py`（含 `__init__.py`），约 3795 行 | **对外接口:** `parent_select.select_parent(info) -> list[dict]`、`parent_select.finalize_parent_ownership(parent, mol) -> dict`、`parent_select.rule_driven_parent_candidates(info) -> list[dict]`、`parent_select.select_principal_parent_skeletons(info) -> PrincipalParentSelection`

---

## 概述

Layer2 接收 [[architecture/layer1-analyzer]] `analyze()` 的 info 字典，按 P-44 选出母体结构 (parent hydride)，返回 **P-45.2.1 并列最优候选组**；每个候选是一个 dict，供 [[architecture/layer3-substituents]] 提取取代基、[[architecture/layer4-numbering]] 在母体骨架上编号、[[architecture/layer5-name-assembly]] 组装名称。

| 文件 | 职责 |
|---|---|
| `parent_select.py` | 层唯一门面：P-44 编排、候选收集、原子归属、P-45.2 排序 |
| `principal.py` | P-41 注册表与主基团选择 |
| `parent_skeleton.py` | P-44 骨架枚举与筛选谓词（含杂原子烃、阳离子/氮负离子短路） |
| `principal_expression.py` | 主基团表达 facts 与母体 dict 合成 |
| `ring_scaffold.py` | 骨架规格 + 保留 SMILES 模板表 + 环解析（唯一事实来源） |
| `kind_registry.py` | kind → `KindMeta`（中英词干/环型/环数/保留名）+ 出口补词干 |
| `hantzsch_widman.py` | P-22.2.2 生成式杂单环：未命中模板的 3–10 元杂环词干/位次/名称与稠合组分名 |
| `fused_system.py` | P-25.3.2.4 稠环拆解为 `FusedNode` 树 |
| `bridged_system.py` | P-23 扩展 von Baeyer 桥环拆解为 `BridgedNode`（主环/主桥/次级桥 + 位次） |
| `spiro_system.py` | P-24 螺环拆解：单环组分走螺描述符（`SpiroNode`），含多环组分走组分式命名（`FbsNode`） |
| `chain_walk.py` | 链行走（等长最长链全枚举，碳与杂原子烃共用同一套元素参数化） |

**输入 (info)**：`mol`（RDKit Mol）、`carbon_ids`/`n_carbons`、`fg_inventory`、`double_bonds`/`triple_bonds`、`rings`/`n_rings`/`has_ring`、`ring_systems`/`n_ring_systems`，以及 `namer._name_mol` 注入的 `root_ctx` 与 `salt`。每个环系 dict 含 `atom_ids`/`sssr_indices`/`fusion_edges`/`spiro_edges`/`free_spiro_atoms`/`n_rings`/`topology`（螺环系为 `"spiro"`）：L1 把螺环对并入同一连通分量，稠合与螺环共用同一个环系。FG 存在性一律由 `fg_inventory` 内容判定；`OXO_FG_CLASSES`、`OXO_CENTER_KINDS`、`FunctionalGroupClass` 全部类别（含 `HETERANE`/`AZANIDE`/`THIOL`）与母体优先序 `PARENT_SENIOR_ATOMS` 由 L1/`constants` 定义，本层只作消费方。

**输出 (parent dict)**：`kind`（母体类型）、`chain`（骨架原子序）、`n_carbons`、`owned_atoms`（frozenset）、`covered_principal_ids`、`principal_occurrences`、`principal_group_count`、`principal_expression_facts`、`stem_en`/`stem_zh`、`numbering_scaffold`、`scaffold_id`/`scaffold_identity`/`scaffold_match`、`hydro_atoms`，环骨架身份字段 `fused_tree` 或 `spiro_node`/`spiro_nodes` 或 `fbs_node`/`fbs_nodes` 或 `bridged_node`/`bridged_nodes`（四者互斥），以及按 FG 与拓扑分支写入的锚点字段、`double_bond(s)`/`triple_bond(s)`、`free_valence_order`、`radical_anchor_element`、`single_atom_skeleton`、`anion`、`o_idx`/`alkoxy_n`/`thio_side`、`hal_idx`/`hal_z`、`oxo_kind`/`n_oh`/`n_om`/`salt_meta`、`amide_z`/`hydrazide_n_idx`/`hydrazide_near_n_idx`，以及杂原子烃与氮负离子分支的 `heterane_z`/`lambda_n`/`lambda_atoms`/`stem_bare_en`/`stem_bare_zh`。

**职责边界**：本层不切割取代基（L3）、不编号（L4）、不拼接名称字符串（L5），也不做烷烃兜底——无候选时返回空列表，由 `namer._run_candidates` 以 `no_assemblable_candidate` 显式失败。反向依赖只用于「借算」：编号候选交 L4 `narrow`、各环系统的编号枚举与稠合/桥环主体名分别向 L4/L5 查询，本层不复刻这些判据。

层的调用链：

```
select_parent(info)
├─ _collect_candidates          rule_driven_parent_candidates → 剔除 None
│   ├─ select_principal_parent_skeletons   主基团选择(P-41) + 骨架枚举与窄化(P-44)
│   └─ _express_selected        express_ring_principal / express_chain_principal
│       └─ _ring_scaffold_and_nodes
│             P-24 螺环短路 → 否则 scaffold 身份 + P-25 稠环树 + P-23 桥环节点（桥环接管身份）
├─ _finalize_ranked             pack_parent_stem(词干/编号 facts) → finalize_parent_ownership(owned_atoms)
├─ _reorder_p45_2(tied=True)    P-45.2.1 前缀取代基团数最大组
└─ _reorder_oxo_ester_side      缩合磷酸酯并列候选取糖/多元醇侧
```

管线内部全部为无副作用纯函数（`tools.memo` 记忆化除外），可独立单测：L2 的输入只有 info，输出只有候选列表，不写回 info。

---

## 核心逻辑

### 骨架枚举与窄化（P-44）

`ParentSkeleton(topology, atom_ids, covered_principal_ids)` 与 `SkeletonSelection(candidates)`；`SkeletonTopology` 取 `ACYCLIC` / `RING_SYSTEM`。

`enumerate_principal_skeletons` 有三条前置短路，均不进通用链/环枚举：

- **阳离子**：occurrences 非空且**全部**为 `FunctionalGroupClass.CATION` 时，只取 `_cation_candidates`（每个锚点一个单原子 `ACYCLIC` 骨架，覆盖锚点所属 occurrence）再过 `keep_senior_atom`——P-73.7(c) 的「多阳离子中心取优先元素」比 P-44 的拓扑规则更专。
- **氮负离子**：全部为 `AZANIDE` 时同样只取 `_cation_candidates`（P-72.2.2.2(2) 氮负离子自任母体，全部臂退为前缀）。
- **杂原子烃**：全部为 `HETERANE` 时取 `_heterane_candidates + _ring_candidates`，再过 `keep_senior_atom`——P-21 杂原子烃以杂原子自任母体（P-41 类 21–39 高于碳）；环候选一并枚举，同级元素由 P-44.2「环优先于链」裁决，避免把含杂原子的环拆成开链。

其余情况返回 `_chain_candidates + _ring_candidates`：

- **链候选** `_chain_candidates` = `_open_chains` + 按原子集去重。`_open_chains` 以主基团锚点（`_anchors`）中的开环锚点为种子，`_all_chains_through` 全枚举等长最长链（平局交下游裁决），再加 `_pair_chains`（`_chain_through_two` 连两锚点的链）；开环锚点无产出时退化为 `_longest_chain` 单候选（纯烃）。`_demoted_leaf_carbons` 把被压为前缀叶的中心碳列为禁行点（P-61.1.3）。`_chain_coverage` 算链覆盖的 occurrence：胺任一臂在链即覆盖（P-62.2），其余类别要求全部锚点在链内；无锚点的 occurrence 不计。
- **非碳母体氢化物候选** `_heterane_candidates`：按锚点元素分组；单原子组只在该原子非标准价（P-14.1.3）或为 P-44.1.2 标准价单核母体氢化物中心（L1 `is_standard_parent_hydride_center`）时出一个单原子骨架；多原子组用 `_all_chains_through`/`_chain_through_two` 在同元素开链子图上枚举链（banned 收窄到该元素的锚点集），支链与取代基退给 L3。链候选不越出杂原子烃网络。
- **环候选** `_ring_candidates`：对 `info["ring_systems"]` 的每个环系构造一个候选，此阶段**不跑** scaffold 识别（身份延迟到表达阶段）。附着判据 `_ring_attaches`——锚点在环内即为附着；`AMINE`/`ALCOHOL`/`THIOL`/`RADICAL`/`KETONE` 只认直接附着（隔碳 SH 只作 sulfanyl 前缀，P-63.1.5），碳锚定的含氧酸/磺酰胺、以及 kind ∈ `_RING_AS_SUBSTITUENT_KINDS`（`boronic`/`phosphonate`）一律返回 False（中心自任母体，环退为取代基，P-65.3.1/P-68.2.1）；其余类别允许经一个邻居间接附着。
- `chain_walk`：全部函数带元素参数 `z`（默认 6）——`_seed_atoms`/`_element_neighbors`/`_longest_chain`/`_component_leaves`/`_all_chains_through`/`_chain_through_two` 共用同一套同元素开链子图逻辑，`_seed_carbons` 只是 `z=6` 的别名；它给最长链种子降集（开链子图为多碳树时只用开链叶，含环/芳香根或叶集为空则回退全碳）；`_component_leaves` 做一次 DFS 出父表与最深叶；`_all_chains_through` 在锚点为端点时从最深叶走向锚点、锚点在链内时取两臂深度和最大的组件组合，平局臂全枚举；`_chain_through_two` 复用 `_component_leaves` 的父表（开链子图是森林，a→b 路径唯一）。

`select_principal_skeletons` 的筛选链：

1. `keep_max_principal_coverage`——保留覆盖主官能团最多的候选；
2. `keep_p44_1_2`——候选拓扑混合时用 `keep_senior_atom`（按 `_SENIORITY` 判最优先元素），拓扑单一时不动作；
3. 拓扑全为 `ACYCLIC` 取 `keep_p44_3`（`p44_3_key`），否则取 `keep_p44_2`（`p44_2_key`，只留环候选）；
4. 末位 `keep_p44_4_unsaturation`（`p44_4_unsaturation_key`，经 `memo.by_key` 记忆）。

各谓词统一走 L4 `numbering_engine.narrow`（`reverse=True` 取键最大者）。

### P-44 谓词速查

| 谓词 / 键 | 条款 | 判据 |
|---|---|---|
| `keep_max_principal_coverage` | P-44 前置 | 保留覆盖主官能团最多的候选 |
| `keep_p44_1_2` / `keep_senior_atom` | P-44.1.1 / .1.2 | 候选拓扑混合时取含最优先元素者（`_SENIOR_ATOMS = constants.PARENT_SENIOR_ATOMS`，序 N>P>As>Sb>Bi>Si>Ge>Sn>Pb>B>Al>Ga>In>Tl>O>S>Se>Te>C；`_SENIORITY` 给序数，L1 的标准价单核判定与它同源） |
| `keep_p44_2` / `p44_2_key` | P-44.2 | 环候选：有杂原子 > N 计数（有 N 时 senior 位取 0）> senior 元素 > 完整环数 > 原子数 > 杂原子数 |
| `keep_p44_3` / `p44_3_key` | P-44.3 | 开链候选：杂原子数 > 原子数（第三位恒 `None`，不参与比较） |
| `keep_p44_4_unsaturation` | P-44.4 | 多重键数 > 双键数（芳香键按半计并入两者） |

`_senior_atom` 取骨架中实际存在的最优先元素；`_ring_count` 数完全落在骨架原子集内的 SSSR 环数。

### 主基团选择（P-41）

`PRINCIPAL_REGISTRY` 由 `layer1.fg_registry.FG_SPECS` 中 `p41 != 0` 的条目派生，每项为 `PrincipalFeatureSpec(PrincipalPriority(p41_class, p43_path), expression, anchor_fields)`。`_spec_from_fg` 只搬运注册表字段，`feature_spec` 为唯一查表入口；`PrincipalExpression` 只用到 `SUFFIX`。

`select_principal_group(inventory, registry=PRINCIPAL_REGISTRY, mol=None)`：取非 `demoted` 条目对应的类别集合（`demoted` 叶 carboxy/cyano 不作主基团候选），过滤掉无规格者，按 `_effective_priority` 取 `min`；选中后取 `inventory.occurrences(选中类)`，得到 `PrincipalGroupSelection(group_class, occurrences)`；无候选时返回 `FG.NONE` 的空选择。`mol` 供阴离子态判定使用（`parent_select` 传入 `info["mol"]`）。

`_effective_priority(group_class, spec, inventory, mol)` 给出候选类的实际 P-41 优先级，三个特例：

- **氧/硫负离子升类**：`FG.ALCOHOL`/`FG.THIOL` 的**全部** occurrence 都 `anion_os`（中心杂原子带负电荷）且分子内无阴离子型酸（`_has_charged_acid`：羧酸根的 `surr_idx` 带负电，或含氧酸 `n_om > 0`）时，升到 `PrincipalPriority(_ANION_OS_P41=4, ...)`——酚盐/醇盐/硫醇盐作特征基团时高于酸（表 4.1 类 4 > 类 7）；中性酸在场不阻断（降为 carboxy/sulfo 前缀）。
- **阳离子让位**：`FG.CATION` 且 `inventory.has_anion` 为真时返回 `PrincipalPriority(_CATION_ANION_GATED=99, ...)`，置底，从而不参与母体竞争（P-41 表 4.1 类 4 阴离子 > 类 6 阳离子）。
- **含氧酸升类**：`FG.OXOACID` 的全部 occurrence 都满足 `oxoacid_is_acid(payload)`（酸式）时升到 `OXO_ACID_P41`，从而排在酯（9）与酰胺（11）之前；否则按注册表登记值参与比较。

它是 `select_principal_group` 的比较键，也是「主基团互斥」的实现——只选单个最高优先级类别。

注册表覆盖的 P-41 等级（`p41` 值）与主基团类的对应关系：自由基与酰基 1、氮负离子（`azanide`）4、阳离子 6（有阴离子时置底）、羧酸 7、含氧酸 9（酸式时升为 `OXO_ACID_P41`）、酯 9、酰卤 10、磺酰胺（`(1,)`）与酰胺 11、腈 14、醛 15、酮与硫酮（`(1,)`）16、醇（`(1,)`）与硫醇（`(2,)`）17、胺 19、杂原子烃（`heterane`）36。`p43_path` 用于同级内的先后，`PrincipalPriority` 以 `order=True` 的冻结 dataclass 承载，可直接作 `min` 的比较键。

### 表达与 kind 正交化

两个入口 `express_chain_principal`（仅 `ACYCLIC`）与 `express_ring_principal`（仅 `RING_SYSTEM`），由 `parent_select._express_selected` 按骨架拓扑分派；两者都只用骨架覆盖的 occurrence（`_covered`）决定 kind 与 multiplicity，但把 `selection.occurrences`（全部主基团 occurrence）原样交给 `_parent_dict` 写入 `principal_occurrences`，供原子归属与 L4 判位次。

`express_ring_principal` 先调一次 `_ring_scaffold_and_nodes`，把 scaffold 身份与已算好的 P-24/P-25/P-23 节点一次交给 `_ring_kind` 与 `_scaffold_fields` 复用；`_ring_kind` 返回 None 时该候选被丢弃（环自由基在无 scaffold 时即此例）。

`_parent_dict` 合成公共字段：`kind`、`chain`(=骨架原子序)、`n_carbons`、`covered_principal_ids`、`principal_occurrences`、`principal_group_count`、`principal_expression_facts`，与各类专属 fields 合并。

`PrincipalExpressionFacts`：`group_class`、`multiplicity`、`relation`（`IN_SKELETON`/`EXOCYCLIC`，按特征原子是否落在骨架内判定）、`occurrence_ids`、`characteristic_atoms`、`anchor_atoms`（官能团原锚点）、`attachment_atoms`（`_skeletal_attachments`：锚点全在骨架外时取其在骨架内的邻居，两者皆空则回退锚点集）。

`_semantic_anchor_fields` 专管 P-14.4(a) 的固定 locant 1 锚点：仅当类别 ∈ `_SEMANTIC_ANCHOR_FGS = {RADICAL, ACYL}` 且锚点唯一时写一个锚点字段，键取 `_anchor_fields(group_class)[0]`；`_anchor_fields` 对 `NONE` 给 `("none_c_idx", "none_c_idxs")`，其余类别取注册表 `anchor_fields`。链骨架与环骨架都把它并进 fields 的最前面。

#### kind 决定

`_chain_kind(group_class, count, occurrences)`：

| 类别 | kind |
|---|---|
| `NONE`（纯烃） | `"alkane"`（仅 count == 0） |
| `ACYL` | `"acyl"`（羰基头为 locant 1，P-65.1.7.2） |
| `RADICAL` | `"radical"` |
| L1 `functional_group_inventory.OXO_FG_CLASSES`（`OXOACID`/`SULFONAMIDE`） | `_oxo_kind_of(occurrences)`：**取自 L1 payload 的 `oxo_kind`**，全部 occurrence 须同一 `oxo_kind`，否则返回 None 表示不支持 |
| 其余 | `group_class.value`（count ≥ 1），故 `CATION` 得 `"cation"`、`AZANIDE` 得 `"azanide"`、`HETERANE` 得 `"heterane"` |

`_ring_kind` 复用同一张表并向 `_chain_kind` 透传 `occurrences`：自由基在无 scaffold 时返回 None（未知杂环无 -yl 词干，显式失败），有 scaffold 时给 `"radical"`；有 scaffold 且类别在 `_FG_CLASSES`（全部注册 FG 类减去 `NONE`）中时取 `_chain_kind`（None 则继续）；环上外环 `-CHO` 单独放行 `"aldehyde"`（P-66.6.1.1.3）；否则 `_resolved_ring_kind`——`scaffold.id` 不在 `carbocycle`/`benzene`/`fused_hetero` 时 kind 即 id（保留模板给 `naphthalene`/`indole`，桥环给 `"bridged"`，两螺环给 `"mono_spiro"`/`"fused_bridged_spiro"`，生成式杂单环给 `"hw_mono"`），这三者或 scaffold 为 None 时一律给 `"alkane"`（芳香多环与纯碳环的 kind 都收敛于此，身份分别由 `scaffold_id`/`fused_tree` 承载）。

#### 阴离子与酯/硫代酯字段

- `_ANION_FLAG_FGS = {ACID, OXOACID, ALCOHOL, THIOL}`：该类别的全部 occurrence 都是阴离子时写 `anion: True`（`_expression_flags` 供链与环共用），L5 据此转 `-ate`/`-olate`。判据 `_is_anion_occurrence` 看 `surr_idx` 周边原子是否带负形式电荷；醇/硫醇的负电荷在中心杂原子本身（`center_idx`），故另并该原子。
- `_ester_o_idx(mol, payload)`：**羰基碳上另连烃基的单键 O 优先，无 O 时取 S**（硫代酯）。
- `_thio_fields(mol, o_idx)`：o_idx 为 S 时写 `{"thio_side": True}`（P-65.6.3.3.7.1，供 L5 换用 thioate 词尾与斜体 S 位次）。`ester_fields`（环骨架）与 `_chain_ester_fields`（链骨架）都补该标记；单 occurrence 时附 `o_idx` 与 `alkoxy_n: 0`，多 occurrence 只写 `o_idx`。
- `_chain_acyl_halide_fields`：单一 occurrence 时从 `surr_idx` 取卤素，写 `hal_idx`/`hal_z`；环骨架的外环酰卤走同一函数。
- `_amide_fields(info, occurrences, fields)`（kind `"amide"` 的链骨架与环骨架共用）：羰基碳上另有 `=N`/`=S` 时写 `amide_z`（`N`/`S`，P-43 类 16/17 硫代酰胺/亚氨酰胺词尾）；否则 `_hydrazide_ns` 判 `C(=O)-NH-N`（`_is_hydrazide_far_n` 认简单胺型远端 N）时写 `hydrazide_near_n_idx`/`hydrazide_n_idx`（P-66.3.1.1 酰肼尾）。

#### 含氧酸字段与盐门控

`_chain_oxoacid_fields(info, occurrences, fields)`（单一 occurrence 才放行）写 `oxo_kind`/`n_oh`/`n_om`：

- kind **不在** `OXO_CENTER_KINDS`（`phosphate`/`phosphonate`/`sulfate`/`boronic`，即中心自任母体者之外的碳锚定磺酸/膦酸等）时直接返回字段、**不做盐门控**。
- kind 为中心自任母体者才校验金属盐：`n_om > 0` 时要求 `salt["n_metal"] == n_om`；`n_om == 0` 时不允许存在金属。通过后附 `salt_meta`。
- 门控不通过返回 None，该候选被丢弃。

`express_chain_principal` 的含氧酸分派条件是 `selection.group_class in OXO_FG_CLASSES`（L1 `functional_group_inventory` 的表，含氧酸与磺酰胺共用）。

#### 不饱和度与自由基/阳离子字段

`_chain_unsat_fields` 把骨架内 C=C/C≡C 写成 `double_bond`/`triple_bond`（单个）或 `double_bonds`/`triple_bonds`（多个）；`stem_bare_en` 在场（杂原子链）时再并入 `_hetero_chain_polys` 的链内杂原子多重键（P-21.2.2，不在 info 的 C=C/C≡C 表内）；`scaffold_id` 为 `carbocycle`、`bridged` 或 `mono_spiro` 时由 `_kekule_ring_dbs` 补回被芳香感知剔除的环内 C=C，`fused_bridged_spiro` 直接跳过（不饱和已含在各组分名内），保留 mancude 母体覆盖的原子集、或有 `fused_tree`、或为生成式 HW 单环（`is_hw_scaffold`）时（`_implied_ring_atoms`）剔除其隐含多重键。

自由基分支：`_mononuclear_radical` 把杂原子锚点收敛为表 2.1 单核氢化物骨架（`atom_ids` 换成单锚点）——锚点带 +1 电荷且元素在 `CATION_FREE_STEMS` 中时取阳离子词干且不并入氧化态/自由价键级，否则按元素把 S 的氧化态（`SULFUR_STEM_BY_OXO`）、N 的自由价键级（`NITROGEN_STEM_BY_FREE_DOUBLE`）、P 的氧化态（`PHOSPHORUS_STEM_BY_OXO`，无对应词干则失败）并入词干，并写 `radical_anchor_element` 与 `stem_en`/`stem_zh`。元素不在单核表内（Si/Ge/Sn/Pb、B 族、卤素…）时回退到 `PARENT_HYDRIDE_STEMS` 的杂原子烃名并写 `heterane_z`/`lambda_n`，链骨架随后把 kind 转 `"heterane"`。碳锚点走 `_radical_free_order`（`_anchor_free_order` 读锚点与 `*` 哑原子之间的键级）：> 1 时写 `free_valence_order`，链引擎据此出 `-ylidene`/`-ylidyne`。链骨架与环骨架共用这套字段。

阳离子分支：`kind == "cation"` 时 `_mononuclear_cation` 取 `cation_parent_names(Z)`，收敛为单原子骨架并写 `stem_en`/`stem_zh`/`single_atom_skeleton`；锚点带 `=O` 时词干前置 `oxo`（oxophosphanium/oxoazanium），否则氧被整段丢弃；名字缺失则候选失败。

氮负离子分支：`kind == "azanide"` 时 `_mononuclear_azanide` 要求单锚点且为 N⁻（形式电荷 −1），收敛为单原子骨架并写 `stem_en="azanide"`/`stem_zh="氮化物"`/`single_atom_skeleton`（P-72.2.2.2(2)）。

杂原子烃分支：`kind == "heterane"` 时 `_heterane_parent` 要求链条元素全同且该元素在 `PARENT_HYDRIDE_STEMS` 内——单核用氢化物名（键数非标准时附 `lambda_n`），多核用 `hydride_chain_stem` 的裸词干（`stem_bare_en`/`stem_bare_zh`，供链引擎拼 `ene`/`yne`）与各原子的 `lambda_atoms`；两者都写 `heterane_z`（P-21.1.2/P-21.2.2）。多核链（`stem_bare_en` 在场）的 `principal_expression_facts.multiplicity` 归 1（整条链是一个母体氢化物，非 n 个主基团）。

### 生成式杂单环（Hantzsch-Widman，P-22.2.2）

`hantzsch_widman.py` 为未命中保留模板的 **3–10 元孤立含杂单环**生成词干与名称，是这一族环的唯一来源，与 `ring_scaffold` 模板表互补。

**身份与接线**：`identity(info, skeleton)` 在原子数 3–10、环内含杂原子、元素都在 `HW_PREFIX_EN` 词表内、且 `isolated_ring`（不与他环稠合、且该环是 SSSR 一员，稠合/桥环另走 P-25/P-23）时返回 `ScaffoldIdentity(HW_ID="hw_mono", HW_CLASS="heterocycle", 1, "hetero")`。`resolve_ring_scaffold` 在模板匹配失败后调用它，故 kind 为 `"hw_mono"`、`scaffold_id`/`scaffold_identity.id` 同为 `hw_mono`。`_scaffold_fields` 对 `is_hw_scaffold` 的骨架改用 `hantzsch_widman.hydro_atoms` 补 `hydro_atoms`（mancude 参照之外的氢位，P-54.4.1）；`_implied_ring_atoms` 把整个环原子集视为母体名隐含不饱和（不再摊成 ene/yne 位次）；`kind_registry.pack_parent_stem` 在 `parent_names(scaffold_id)` 未命中后调 `hw_parent_names` 取词干。

**词干**：`sat_stem`/`unsat_stem` 按环大小与 `six_group` 出 Table 2.5 词干——三元 `irane/irine`（全氮 `iridine`）、四元 `etane/ete`、五元 `olane/ole`（含氮 `olidine`）、六元 `HW_SAT_SIX`/`HW_UNSAT_SIX`（A 组 `ane`、B 组 `inane`、C 组 `inine`）、七至十元 `HW_SAT_TAIL`/`HW_UNSAT_TAIL`。`six_group` 由优先性最低的杂原子所属组决定。中文 `zh_stem`：含氮五/六元环取「唑/唑烷」「嗪/嗪烷」（N 折入基干，`_zh_folded_n` 加倍数词），其余按 `HW_ZH_RING` 环前缀 + `_zh_unsat_tail`。

**位次与省略**：`locant_map` 按 P145_SENIOR 引用顺序分组杂原子位次，`locant_string` 按该顺序（非数值升序）成串；`omit_locants` 依 P-22.2.2.1.7——单杂原子、或 `_arrangement_classes`（杂原子多重集在环上的循环排布数，模旋转翻转）为 1 时省略全部位次。`elide_a` 按 P-22.2.2.1.1 做前项尾 `a` 的元音省略（`tetra`+`aza`→`tetraza`）。`ring_lambda_atoms` 对 `HW_RING_LAMBDA_Z`（P/As/Sb/Bi/B/I/Si/Ge/Sn/Pb）中键数偏离标准值的杂原子在位次后附 λ（P-22.2.7.1）。`ring_numbering` 复用 L4 `_narrow_hetero_ring` 编号，保证名中位次与 L4 同源。`hw_name_from_cycle(zs, n_db, lambda_at)` 是组装核心，`parent_names(sid, mol, chain)` 是其对外入口。

**加氢**：`mancude_hydrogens` 用 `_max_matching_mask`（k≤10 直接枚举边子集）求环内最大匹配得参照氢数与双键数，`hydro_atoms` 取实际氢数多于参照者；环内无重键时取饱和词干、不产生加氢前缀。

**稠合组分**：`component_key(mol, atoms)` 给 3–10 元、含杂、mancude（环内双键数等于参照）的孤立环编出 `"hw:" + 环序元素符号串`（`HW_COMPONENT_PREFIX`，P-25.2.2.1.1）；`component_names(sid)` 反解为稠合组分词干，位次非空时加方括号（P-25.3.2.1.2）。`ring_scaffold.match_fusion_component` 在保留母体与单环烃之后接 `component_key`，`component_stem` 以 `hw:` 前缀分派到 `component_names`，故未注册杂环同样能作稠合零件。

### kind 与词干注册

`kind_registry` 是母体 kind 的中英词干权威，注册表 `_REG` 只由 `ring_scaffold.all_specs()` 经 `_load_from_scaffold_specs()` 在 import 期填充：每个带 `stem_en`/`stem_zh` 的 `ScaffoldSpec` 生成一条 `KindMeta(kind, en, zh, ring, n_rings, retained)`。`get`/`parent_names` 是唯一读入口；未注册 kind 返回 None。`mono_spiro`/`fused_bridged_spiro`/`bridged`/`hw_mono` 等非模板身份不在注册表内：桥环/螺环词干由 L5 的组装器现算，`hw_mono` 词干由 `hantzsch_widman.parent_names` 现算。

`pack_parent_stem(parent, mol)` 在 L2 出口补齐词干与编号字段：

- 按 `parent_names(scaffold_id)` → `hw_parent_names(scaffold_id, mol, chain)` → `parent_names(kind)` 三级回退取名；命中且尚未写 `stem_en`/`stem_zh` 时写入。
- 再由 `locant_prefix(scaffold_id)` 决定是否前置 locant 前缀（如 `1H-`、`1,3-`）；词干自带同一前缀时不重复前置；条件化前缀（`prefix_nh_conditional`）只在该环有未取代芳香 NH 时保留（`_ring_keeps_nh_prefix`，按 `chain` 中芳香 N 的 H 数判定），否则去前缀。
- 末步 `_attach_numbering_scaffold` 写入 `numbering_scaffold`（固定编号标签表），解析失败则原样返回。
- `mol` 缺省时把当前 mol 合并进候选，供词干判定取原子状态。

### 环骨架模板表（`ring_scaffold`）

`ring_scaffold.py` 是骨架身份的唯一事实来源。

**数据结构**：`NumberingPolicy(standard_path)`、`ScaffoldSpec(id, naming_class, stem_en, stem_zh, n_rings, ring, retained, numbering, locant_prefix, prefix_nh_conditional)` 与 `ScaffoldIdentity(id, naming_class, n_rings, ring)`；`ScaffoldSpec.identity` 给出身份。

**模板表**：`_TEMPLATES` 共 **103 条**保留 SMILES 模板，每条含 `smiles`/`naming_class`/`stem_en`/`stem_zh`，另有可选 `fused`、`fused_prefix`、`fused_stem`、`locant_prefix`、`prefix_nh_conditional`、`standard`。命名类分布：`monohetero` 50、`naph_family` 19、`fused56` 12、`phenothiazine` 3、`purine` 2、`xanthene` 2，以及 `mono_carbo`/`anthra`/`phenanthrene`/`pyrene`/`chrysene`/`picene`/`carbazole`/`acridine`/`benzodioxole`/`indolizine`/`pyrrolizine`/`steroid`/`pentalene`/`phenalene`/`adamantane` 各 1。`fused` 为假的模板（`thiopyran`、`dithiole`、`dithiolane12`、`adamantane`、各饱和单杂环）只作游离母体，不作稠合零件；`phenalene` 无 `fused` 同样不作零件。未登记进模板的 3–10 元饱和/不饱和杂单环由生成式 HW 路径命名（见上），模板只留 IUPAC 保留名（如 `pyrazolidine`、`imidazolidine`、`dithiolane12`、`thiomorpholine`）。

`_spec_from_template` 由模板派生 `ScaffoldSpec`（环数与环型由查询分子自动算、`retained=True`、`numbering.standard_path` 取 `standard` 标签），`all_specs()`/`get_spec()` 对外供全库使用；`kind_registry._load_from_scaffold_specs()` 在 import 期据此建 `KindMeta` 词干权威，`pack_parent_stem` 按 `scaffold_id` → 生成式 HW → `kind` 三级回退取名。

**固定编号**：49 条模板携带 `standard = (labels, order)`，含 `dihydrofuran`/`dihydropyran`/`dihydropyrrole`/`dihydroimidazole`/`dihydrothiazole` 这类部分不饱和环（字面位次即固定编号，P-14.4(a)/(b)）与传统编号骨架（蒽/菲/芘/咔唑/吖啶/吩噻嗪/噻蒽/呫吨/甾体/金刚烷等）。import 期由 `_STANDARD_LABELS`/`_STANDARD_ORDERS` 两个推导表承载；`standard_chain` 把固定编号映射为分子原子序，`numbering_scaffold_facts` 在标签数与骨架原子数一致时物化 `numbering_scaffold`。

**匹配与身份**：`match_retained`/`_match_with_map` 先按元素签名（`_elem_sig`/`_TEMPLATE_ELEM`）剪枝，再要求子图同构原子集**精确等于**骨架原子集，并由 `_is_induced_match` 校验原子集诱导子图的键数等于模板键数——子图同构容忍目标多出的键，缺此校验则环系的真子图模板会把多出的环静默丢掉。`_isolated_saturated_ring` 以 `hantzsch_widman.isolated_ring` 判孤立单环、再判 Kekulé 视图环内是否全为单键，是则只允许匹配 `_HAS_MULTI_BOND` 为假的模板（P-31.2）；`mancude_only` 时跳过 `fused` 为假的饱和保留名（不作稠合零件）。`_match_with_map` 内部用 `_scan` 先在原分子上扫 `_Q`，未命中再在完全氢化的分子副本（`_hydrogenated`，经 `memo.by_mol` 记忆）上扫 `_Q_H`（P-25.3.4），两轮共用同一套元素签名与诱导覆盖校验；氢化轮命中 `_HAS_MULTI_BOND` 为假的模板时，若分子在该原子集内仍含非单键（按饱和保留名会整段丢掉不饱和度），返回 None 改走生成式 HW。`resolve_ring_scaffold` 命中则返回 `spec.identity`；否则交给 `hantzsch_widman.identity`（未命中模板的 3–10 元杂单环）；再无则 `_generic_carbocycle`——按骨架内完整环数（`parent_skeleton._ring_count`）返回 `fused_hetero`（≥ 2）/ `carbocycle`（单环，`ring` 记 `carbo`）。`_spec_from_template` 的环型由 `_ring_kind` 统一判定。

**命名类速查**：

| naming_class | 条数 | 代表模板 | kind 归属 |
|---|---|---|---|
| `monohetero` | 50 | 呋喃/吡啶/吡喃/噻喃/噁唑/哌啶/硫代吗啉/氧杂环己烷 | kind 即 scaffold id |
| `naph_family` | 19 | 萘、喹啉、异喹啉、色烯、蝶啶、酞嗪、苯并二噁英、苯并二氮杂卓 | kind 即 id |
| `fused56` | 12 | 吲哚、苯并呋喃（1-/2-）、苯并噻唑、苯并三唑、苯并二氧杂环戊烯 | kind 即 id |
| `xanthene` | 2 | 呫吨、噻吨 | kind 即 id |
| `phenothiazine` | 3 | 吩噻嗪、吩噁嗪、噻蒽 | kind 即 id |
| `purine` | 2 | 嘌呤、吡唑并[5,4-d]嘧啶 | kind 即 id |
| `mono_carbo` | 1 | 苯 | kind 收敛 `alkane` |
| `carbocycle`（非模板） | — | 未注册单环（无模板命中且非含杂时按环数兜底） | kind 收敛 `alkane` |
| `fused_hetero`（非模板） | — | 未注册稠环 | kind 收敛 `alkane`，身份由 `fused_tree` 承载 |
| `bridged`（非模板） | — | 已注册或未注册的桥环 | kind `bridged`，身份由 `bridged_node` 承载 |
| `mono_spiro` / `fused_bridged_spiro`（非模板） | — | P-24 螺环系 | kind 同名，身份由 `spiro_node`/`fbs_node` 承载 |
| `hw_mono`（非模板） | — | 生成式杂单环（P-22.2.2） | kind `hw_mono`，词干由 `hantzsch_widman` 现算 |
| 其他单例 | 各 1 | 蒽/菲/芘/苉/屈/咔唑/吖啶/甾体/戊搭烯/菲那烯/金刚烷/中氮茚/吡咯嗪 | kind 即 id |

**加氢与指示氢**：`hydrogenated_atoms`（模板 Kekulé 双键位在分子中已饱和者；环杂原子与碳位的新增 H 同为加氢位，一并进 hydro 前缀，P-31.2.2；计数落在 `HYDRO_MULT_N` 之外且存在环外多重键时，其中一位改由指示氢承载，P-58.2.1）、`extra_indicated_atoms`（模板同位无 H 而分子有 H 的芳香杂环原子，仅稠合母体）、`extra_hydrogenated_atoms`（氢数多于模板同位的环原子，或模板双键位已饱和且 H 数不增——该位另有取代基，如 4H-异喹啉-1,3-二酮的 C4）、`mancude_ring_atoms`（保留 mancude 母体覆盖的整个不饱和环）、`_kekule_double_atoms`（模板双键端点缓存）。`locant_prefix` 返回 `(en, zh, nh_conditional)`。

**稠合前缀**：`retained_fusion_prefix` 先查模板 `fused_prefix`，否则回退 `fusion_carbocycle_prefix`；`component_stem` 给出组分的稠合词干（`fused_stem` 优先，否则经 `_bracketed_stem`）。要点：

- 饱和环作稠合组分须取 **mancude 对应名**（P-25.3.1.2.1 / .2.3），而不是把饱和环名去尾加「并」：`oxane → ("pyrano", "吡喃并")`、`thiolane → ("thieno", "噻吩并")`、`thiane → ("thiopyrano", "噻喃并")`；`azepane` 的 `fused_stem` 为 `("azepine", "氮杂卓")`（P-25.3.1.2.2）。
- `_bracketed_stem`：杂原子位次前缀为纯数字（`1-`/`1,3-`/`1,2,4-`）且词干本身不以数字开头时，稠合组分名取方括号形式（`[1,3]oxazole`、`[1]benzofuran`、`[1,3]benzodioxole`），即 P-25.3.5 的「稠合母体须引用完整位次」；`1H-` 型前缀与词干已带位次者（`1,2-oxazole`）不加方括号。
- `_FUSION_CARBOCYCLES`（6 条单环烃附加组分，P-25.3.2.2.1）：`cyclohexane` 的稠合前缀是保留前缀 **`benzo`/`苯并`**（该条「除 benzo 外」的例外），其余为 `cyclopropa`/`cyclobuta`/`cyclopenta`/`cyclohepta`/`cycloocta`。查询按 SMARTS 现造，键级用 `~` 通配——附加组分只论环大小与元素，不饱和碳环（环戊烯/环己烯等）同样可作稠合零件。`match_fusion_carbocycle` 做精确匹配，`omits_fusion_numbers` 判 benzene 与单环烃附加组分省略数字位次（P-25.3.8.1）。

### 稠环树（P-25）

`_ring_scaffold_and_nodes` 取与骨架原子集一致的环系 dict（`_skeleton_system`），其 `sssr_indices` ≥ 2 时先调 `decompose_fused_system(info, system)` 得 `FusedNode` 树，再调 `try_bridged_scaffold` 得 P-23 桥环节点。

`FusedNode`：`scaffold_id`、`atom_ids`、`ring_indices`、`fusion_shared`（与母体的稠合共享原子）、`attached`（附加组分，递归）、`fused_stem`、`fused_prefix`、`fused_omit_numbers`。

拆解流程：`_candidates_for` 从每个自身可作稠合零件的环种子（`match_fusion_component` 命中）出发，沿融合图（`_fusion_adj`）按栈增长枚举候选环集，按原子集去重。`_select_base` 按 P-25.3.2.4 准则 (a)–(j) 逐条 `narrow` 收窄——(a) 最优先杂原子、(b) 环数多、(c) 环大小降序、(d) 杂原子总数、(e) 杂原子种类、(f) 最优先杂原子计数、(g) 水平行环数、(h) 杂原子位次最低，逐元素位次最低，稠合碳位次最低；准则 (g)–(j) 依赖 L4（`preferred_orientations` + `number_fused_system`，经 `_numbered_locants`）算位次，可编号候选不足两个时该步原样返回。全部并列时取环集升序最小。`_decompose` 递归：选出母体组分后，剩余环按融合图连通分量（`_ring_components`）递归为附加组分，共享原子记入 `fusion_shared`。

`decompose_fused_system` 的两个前置否决：环系带 `free_spiro_atoms` 时返回 None（螺连结的环组分不参与稠合拆解，留给 P-24）；整环系本身已是一个 `component_stem` 为 None 的保留母体（如金刚烷这类笼状桥烃）时返回 None——稠合拆解对它无意义且会误判加氢。

`match_fusion_component` = `match_retained(mancude_only=True)` 或 `match_fusion_carbocycle` 或 `hantzsch_widman.component_key`（生成式 `hw:` 组分），按 `(id(mol), 环原子集)` 经 `memo.by_key` 记忆（对称笼架多个环集张成同一原子集）。

### 桥环节点（P-23）

`bridged_system.py` 表达 P-23 扩展 von Baeyer 桥环：把环系拆成主环、主桥与次级桥，并给出全原子位次。

**数据结构**：`BridgeSegment(heads, atoms, numbers, locants)` 描述一段桥，`heads` 为两端桥头、`atoms` 沿 `heads[0] → heads[1]` 序；`__len__` 返回桥内原子数，即 von Baeyer 描述符中的数字（P-23.2.6.1.2），`reversed()` 掉转方向。`BridgedNode(atom_ids, numbering, main_ring, ring_segments, main_bridge, independent_bridges, dependent_bridges, ring, n_rings)` 是桥环母体：`main_ring` 为编号序的主环原子、`ring_segments` 为两条环段的长度、`ring` 取 `carbo`/`hetero`。只读派生属性：`secondary_bridges`（独立桥在前、依赖桥在后，即引用顺序）、`descriptor`（主环两段长 + 主桥长 + 各次级桥长）、`locant_pairs`（次级桥的上标位次对，小者在前，按引用顺序）；`scaffold_identity()` 返回 `ScaffoldIdentity("bridged", "bridged", n_rings, ring)`，id 与命名类同名。

**识别流程**：`decompose_bridged_system(info, system)` 是公共入口，非桥环返回空表。前置筛：原子数 ≥ 4、`_ring_count`（键数 − 原子数 + 1，P-23.2.6.1.1）算出的环数 ≥ 2、桥头 `heads` 取度 ≥ 3 的原子（P-23.1.1）且至少两个。`_bridges` 拆出全部桥：每个「去桥头分量」须为简单路径且两端各邻接一个不同桥头（`_path_segment`，单原子分量被两桥头夹持时即为 1 原子桥），再加桥头之间的直键（0 原子桥，P-23.1.2）；并要求 `len(bridges) - len(heads) + 1` 等于环数。`_candidates` 在桥图 H（顶点 = 桥头、边 = 桥，`_h_graph`）上枚举简单环（`_cycles`，上限 `_MAX_CYCLES`，两顶点间的平行边视为二元环）作主环；主桥由 `_main_paths` 枚举——H 中两端落在主环上、内部顶点落在环外的简单路径（P-23.1.2/P-23.2.4），链内原子本身可为桥头（如 phenalene 型中心原子），此时并成一条主桥，故 `_Cand.main_e` 是桥边元组而非单边；每对端点生成正反两个方向。其余桥须构成以主环或主桥原子（`core = rset | 主桥环外顶点`）为端点的简单链（`_secondary_ok`：core 外顶点度恰为 2）。`_to_node` 给候选编号，编号失败返回 None；候选组装后还要求描述符数字之和加 2 等于骨架原子数（P-23.2.6.1.4），不符者丢弃。环型由 `ring_scaffold._ring_kind` 统一判定。

**收窄阶梯** `narrow_candidates` 依次施加（每步走 L4 `narrow`，末步不做确定性 tie-break）：主环原子数最大（P-23.2.1）→ 主桥最长（P-23.2.4）→ 主环两段长之差最小（P-23.2.6.2.1）→ 独立桥长度序（P-23.2.6.2.2）→ 依赖桥最少（P-23.2.6.2.3）→ 上标位次集合最低（P-23.2.6.2.4）→ 引用序位次最低（P-23.2.6.2.5）。

**编号**：主环从长段一端进入——长边在前（P-23.2.3），必要时把两条环段的行进方向一并掉转以保持路径连续；主环依次占位次 1..n，主桥由 `_main_bridge_atoms` 沿 `main_e` 展开桥内原子（链内桥头本身也是主桥原子）自桥头续编，`_number_one` 从紧邻较高编号桥头的一端开始给桥内原子编号（P-23.2.6.3），`locants` 记两端桥头的位次对。次级桥分独立桥（两端桥头都已入编号，即在主环或主桥上）与依赖桥。`_number_secondary` 中独立桥的**编号顺序**按桥头最大位次降序（P-23.2.6.3），但返回值仍保持调用方给的**引用顺序**；依赖桥逐个取两端桥头皆已编号者，出现不可编号者则该候选失败。引用顺序由 P-23.2.6.2.2/.5 裁决：`_orders` 固定长桥在前、只对同长度桥排列，枚举数超过 `_MAX_ORDERS` 时退化为「长度降序」单一顺序；逐个顺序编号后取（次级桥长度降序、引用序位次展平）键最小者。

**路由** `try_bridged_scaffold(info, scaffold, fused_tree, system)`：环系带 `free_spiro_atoms` 时直接让位（归 P-24）；再取 `decompose_bridged_system`，空表即让位。scaffold 为 None 或 id 非 `fused_hetero`（保留模板与环系同构，P-25.2）时维持稠环路径、不走桥环；`_has_retained_peri_parent` 命中（稠环树里含芘/phenalene 这类保留的迫位稠合母体，其环集的环邻接图非树）时同样维持稠环路径（P-25.1.1，保留母体优先）。**条件 A** —— `_fusion_naming_applies` 为假时直接走桥环，该判据要求：环邻接图（顶点 = SSSR 环、边 = 共享成键的两环）为一棵树（`_ring_adjacency_is_tree`，P-25.5/P-52.2.4.4——非树或非连通时「组分对」不再唯一，稠合名原理上不可行）、环系内至少两个 ≥ 5 元环（P-52.2.4.1）、且 `_cannot_be_mancude` 为假（环系内已成满 4 根 σ 键的碳腾不出 π 键，写不出最大非累积双键，P-25.3.1.2）。否则要求 `_tree_covers_rings` 确认稠环树的 `ring_indices` 覆盖环系的全部 `sssr_indices`（盖不住环系的全部环则回退桥环）；末步 `_p25_names_ok` 反问 `layer5.fused_namer.fused_parent_names`——**能拼出稠合词干（条件 B）则保持稠环路径**，拼不出才回退桥环。

桥环命中时以 `bridged[0].scaffold_identity()` 接管 scaffold 身份，`kind` 由 `_resolved_ring_kind` 给出 `"bridged"`。两套拆解都独立于 scaffold 模板：未注册系统同样产出。

### 螺旋环系统（P-24）

`spiro_system.py` 只在环系带**自由螺原子**时接管。自由螺原子由 L1 `ring_systems` 判定：从环子图去掉该原子后断成 ≥2 个分量（P-24.1）；`spiro_edges` 为 `(环 i, 环 j, 螺原子)` 三元组表。`SPIRO_SCAFFOLDS = ("mono_spiro", "fused_bridged_spiro")` 是本层两个螺环身份 id。

#### 与稠环/桥环的互斥路由

`_ring_scaffold_and_nodes` 是三条路线的唯一分岔口，返回四元组 `(scaffold, fused_tree, bridged, spiro)`。分岔顺序固定：**P-24 先于 P-25/P-23 判定**——只要环系的 `free_spiro_atoms` 非空就整体短路到 P-24，`fused_tree` 与桥环节点都不参与计算；否则才做 scaffold 识别并依次求稠环树与桥环节点（后者可接管 scaffold 身份）：

| 环系特征 | scaffold 身份 | 节点字段 |
|---|---|---|
| `free_spiro_atoms` 非空 | `spiro_scaffold_identity`：全为单环组分 → `mono_spiro`，否则 `fused_bridged_spiro` | `fused_tree` 为 None、桥环为空表；`mono_spiro` 写 `spiro_node`/`spiro_nodes`（节点来自 `decompose_spiro_system` 的 `SpiroNode`），`fused_bridged_spiro` 写 `fbs_node`/`fbs_nodes`（节点来自 `decompose_fbs_system` 的 `FbsNode`） |
| 无自由螺原子，`sssr_indices` ≥ 2，稠合命名法适用 + scaffold 身份为 `fused_hetero` + 非保留迫位稠合母体 + 稠环树盖满全环 + P-25 拼得出词干 | `resolve_ring_scaffold` 的结果（保留模板 id 或 `fused_hetero`） | `fused_tree` |
| 无自由螺原子，`sssr_indices` ≥ 2，其余情况（稠合不适用 / 身份非 `fused_hetero` / 保留迫位稠合母体 / 树盖不住 / P-25 出不了词干） | 桥环节点的 `scaffold_identity()`（`bridged`）接管 | `bridged_node` + `bridged_nodes` |
| `sssr_indices` 只有 1 个环 | `resolve_ring_scaffold`（保留模板 id、`hw_mono` 或 `carbocycle`） | 三者皆无 |

两张守卫在 P-24 之外也独立成立：`decompose_fused_system` 与 `try_bridged_scaffold` 各自在 `free_spiro_atoms` 非空时返回空。`_scaffold_fields` 按互斥优先级写字段——螺环身份 → 桥环 → 稠环树，`scaffold=None` 时该函数自行调一次 `_ring_scaffold_and_nodes`。

#### 全单环组分：P-24.2 螺描述符

`_segments(rings, indices, spiros)` 逐环分解出段：环上**无**螺原子 → 返回 None，表示环系含非螺连的稠合环（转组分式命名）；恰 1 个螺原子 → 端环，整环自该原子起绕一周为一段（`terminal=True`，段的 `a == b`）；≥2 个 → 中心环，环序上相邻螺原子之间各成一段。`_cycle_from` 做环序旋转，`_flip` 掉转段方向得到 `_Seg(ring_index, a, b, atoms, terminal)`。

编号候选的枚举与编号是两趟：

- `_enumerate_walks(segs, cap=_MAX_WALKS=64)` 从端环段（无端环则全部段）出发、正向与反向各试一次，`_step` 深度优先串起全部段成回路。每步的候选顺序为「端环段优先（键 0）→ 中心环段（键 1）、再按段内原子数升序」，即 P-24.2.2 的最短路径优先；回路数达上限即停。
- `_number_walk(segs, walk)` 按段引用顺序编号：先给段内原子依次编号并记段长入 `descriptor`，再看段末端 `seg.b`——已编号则把该位次记入 `descriptor_superscripts`（重访螺原子，P-24.2.2 的上标），未编号则新占一个位次、上标记 0。
- `len(free) == 1`（单螺）时把全部上标置 0——P-24.2.1 的单螺描述符不带上标，上标是多螺专用。

`SpiroComponent(index, atom_ids, ring_index, spiro_atoms)` 记每个环组分。候选经 `_numbers_all` 过滤：`numbering` 必须覆盖骨架全部原子且取值恰为 `1..n` 无重号。存活者组装为 `SpiroNode(scaffold_id="mono_spiro", atom_ids, free_spiro_atoms, components, descriptor, descriptor_superscripts, numbering, ring, n_rings)`，`scaffold_identity()` 给 `(mono_spiro, mono_spiro, n_rings, ring)`。全部段回路都是并列候选，本层不裁决，编号选择交 L4 `spiro_numbering`（按 P-24.2.2/2.3 螺位次、描述符数字、杂原子、后缀、双键、取代基次序收窄）。

#### 含多环组分：P-24.5~24.7 组分式命名

`decompose_spiro_system` 是公共入口：`free_spiro_atoms` 为空返回空表；`_segments` 为 None 时转 `decompose_fbs_system`（`FbsNode`），否则出 `SpiroNode` 列表。

- `component_indices(system)` 用并查集把环下标切成组分：`fusion_edges` 一律合并，`spiro_edges` 只在螺原子**非**自由螺原子时合并（非自由螺原子不切断，两侧仍属同一多环组分）；按最小环下标排序返回。
- `_component_system(system, comp, atoms)` 造只含本组分的环系 dict（`atom_ids` 换成组分自身原子，否则环数会算错；`spiro_edges`/`free_spiro_atoms` 清空），组分内部照常走 P-23/P-25。
- `_build_component(mol, info, system, comp, free, k)` 给单个组分定基名与编号候选：
  - **保留名组分**（`match_retained(info, atoms)` 命中）：`kind` 取 `mono_ring`（`len(comp) == 1`）或 `fused_ring`。`_retained_matches` 取模板在组分原子集上的**全部**映射（先 `_Q`，无则退完全氢化模板 `_Q_H`，再无则取 `_match_with_map` 的主映射），对称保留母体的自同构因此各得一个取向候选（如 2-benzofuran 的 1/3 位互换），哪个取向入选由 L4 按 P-24.5.2 的螺位次最小定。`_retained_numberings` 优先用 `_STANDARD_LABELS`/`standard_chain` 的固定编号视图，否则单环走 `_mono_numberings`（L4 `_ring_cands` 全旋转/翻转，杂环再经 `_narrow_hetero_ring`，P-22.2.2.1.3），多环走 L4 `fused_component_numbering`（稠合点作 `sub_layer` 逐层最小化，P-25.3.1.3）。基名由 `_retained_base` 调 `pack_parent_stem` 取（含 `a` 位次前缀，与母体出口同口径）；非单环时把基名里的方括号片段抽成 `order_extra`（P-24.5.3 的斜体稠合字母平局键）；组分内多氢的饱和位由 `_component_hydro_prefix` 补加氢前缀（`hydro` 字段）。
  - **未注册单环**：杂环本层不支持（返回 None）；碳环取 `_cyclo_base` 给的全名与裸词干（`cyclo` + 烷词干，裸词干可接 `a-1,3-diene` 型后缀），编号走 `_mono_numberings`。
  - **未注册多环**：`_flat_numberings` 先试 P-23——`decompose_bridged_system` 在组分子环系上出节点、且 L5 `bridged_body_names(node)` 出主体名与裸词干，则 `kind = bridged_ring`、`order_extra` 取 `descriptor`、编号取各并列节点的 `numbering`（'a' 前缀按 P-24.5.2 提到 spiro 之前，故组分名取裸词干）；否则试 P-25——`decompose_fused_system` + L5 `fused_parent_names` 出双名，编号由 `fused_component_numbering` 给，`order_extra` 取英文名里的方括号片段；两条路都不通、渲染不出基名或编号缺失则返回 None。
  - 基名或编号候选任一缺失 → 组分不可命名 → `decompose_fbs_system` 整体返回空表（显式失败，不静默丢组分）。
- `FbsNumbering(chain, labels)` 是一个组分的**本位次**编号候选（不含撇号），`locants` 属性给「原子 → 位次标签」；`FbsComponent(index, kind, atom_ids, ring_indices, spiro_atoms, base_en, base_zh, bare_en, bare_zh, sid, node, numberings, order_extra, hydro)` 里的 `index` 是引用序槽位（也就是撇号数），`node` 是组分的 `BridgedNode`/`FusedNode`/None，`hydro` 是该组分的加氢前缀 `(en, zh)`（`_component_hydro_prefix`：环内饱和位——螺原子除外——按本组分编号取位次，仅对 `_HAS_MULTI_BOND` 为真的保留名组分、且位次可解析时给出；引用序仍按裸基名算，P-24.5.1）；`FbsNode(atom_ids, free_spiro_atoms, components, links, cites, ring, n_rings)` 是母体，`scaffold_id` 恒 `fused_bridged_spiro`。

引用序 `citation_order(comps, edges)`（`edges[k]` 是组分 k 的 `(螺原子, 对端组分)` 邻接表）：

- 组分数 < 2 或无端组分（度 1）→ None；
- 两组分（P-24.5.1）：`_order_key` 较小者列首；
- 全组分度 ≤ 2（直链，P-24.6）：从 `_order_key` 较小的端组分出发沿链展开，遇到非单一路径则失败；
- 否则要求恰有一个度 ≥ 3 的中心、其余全为端组分（更深的支链归 P-24.7.4，本层不做）：端组分基名全同时中心列首、全部端组分合并为一个引用组（P-24.7.1，即 `tris(...)` 型）；否则 `_terminals_then_center`——`_order_key` 最早的端组分先列、随后中心、其余端组分按字母序（同名合并为一组）（P-24.7.2，即 `bis(...)` 型）。

`_order_key(comp) = (alpha_order_key(base_en), order_extra)`——先字母数字序，再平局键（斜体稠合字母 / von Baeyer 描述符）（P-24.5.3）。

`decompose_fbs_system` 的组装：组分图的边数须等于组分数 − 1（组分图必须是树，否则失败）；`citation_order` 给出引用项分组（每项长度 1 为单组分、> 1 为倍增词组）；展平后的槽位号覆盖组分 `index`（撇号数）；`links` 按列出顺序逐个给出 `(螺原子, 父槽位, 本槽位)`，每个后列组分必须接在已列组分上，接不上则失败；`cites` 记每个引用项覆盖的槽位；`ring` 由全环系原子是否全碳给 `carbo`/`hetero`。

#### 字段与下游

`_scaffold_fields` 按 scaffold 身份选键名，两条并行字段互斥：`scaffold.id == "mono_spiro"`（`decompose_spiro_system` 出的单环组分螺环）写 `spiro_node`（首个并列 `SpiroNode`，空候选写 None——L4 据此显式失败，不得下沉 P-25）与 `spiro_nodes`（全部并列候选）；`scaffold.id == "fused_bridged_spiro"`（`decompose_fbs_system` 出的组分式螺环）写同构的 `fbs_node`/`fbs_nodes`，元素类型换成 `FbsNode`。两者都与 `bridged_node(s)`、`fused_tree` 互斥。`_chain_unsat_fields` 对 `fused_bridged_spiro` 直接返回字段（不饱和已含在各组分名内，母体层不出 ene/yne 位次），对 `mono_spiro` 与 `carbocycle`/`bridged` 一样用 `_kekule_ring_dbs` 补环内 C=C。

编号与词干全部委派：L4 `orient_numbering` 见到 `fbs_nodes` 走 `fbs_numbering`（逐组分收窄 + `_joint_pick` 按引用序位次元组最低选组合，P-24.6.1/P-24.7.2），见到 `spiro_nodes` 走 `spiro_numbering`，两者都直接返回、无论成败均不下落到通用路径（`_fused_numbering` 的 `chain_fused` 判据对螺环同样恒真，会误吞）；L5 `assembler._ensure_parent_stem` 依次看 `fbs_node` → `spiro_node` → `bridged_node`，分别调 `spiro_namer` 的 `fbs_parent_name`/`spiro_parent_names`。L5 不 import 本层，节点按鸭子类型读（`descriptor`/`components`/`links`/`free_spiro_atoms`）；`spiro_system` 反而向 L5 查询组分名渲染与烷词干。

### 原子归属

`finalize_parent_ownership(parent, mol)` 生成不可变 `owned_atoms` frozenset（已是 frozenset 则原样返回），它是 L3 取代基块的边界（见 [[concepts/atom-ownership]]）：

- **`parent["kind"] in OXO_CENTER_KINDS`**（`phosphate`/`phosphonate`/`sulfate`/`boronic`，中心自任母体）：`owned_atoms` **只取主官能团特征原子、不含链**——链仅供编号，臂一律退为取代基。
- 其他 kind：`owned_atoms` = 链原子 ∪ 主官能团特征原子。

`_kind_fg_atoms` 的算法：

1. 种子 = occurrence 锚点 ∩ 链；种子为空（exocyclic，如苯甲酸的羧基）时改取与链相邻的锚点。
2. **`facts.group_class in OXO_FG_CLASSES`（L1 `functional_group_inventory` 的表，= `{OXOACID, SULFONAMIDE}`）** 时，只并「本骨架覆盖的 occurrence」（`covered_principal_ids`）的 `characteristic_atoms` 后返回——**不沿线借臂**，因为中心 P/S 不是碳骨架成员、没有可外借的臂。
3. 其他类别从种子出发，沿邻接并入属于 occurrence 特征原子、且本身不是锚点的邻居（沿线借臂）。

### P-45.2 排序

`select_parent(info)` = `_collect_candidates` → `_finalize_ranked` → `_reorder_p45_2(tied=True)` → `_reorder_oxo_ester_side`：

- `_collect_candidates`：`rule_driven_parent_candidates`（= `select_principal_parent_skeletons` + `_express_selected`）剔除 None。
- `_finalize_ranked`：逐个候选先 `pack_parent_stem(c, mol)`（补 `stem_en`/`stem_zh`、locant 前缀与 `numbering_scaffold`），再 `finalize_parent_ownership` 固化 `owned_atoms`。
- `_reorder_p45_2`：候选不足 2 个直接返回；排序键为三元组 **`(前缀取代基团数, 缩合磷酸链内桥氧数, 原始序)`**，前两项降序、原始序升序稳定排序；`tied=True`（层出口）只保留第一项并列最大的一组，即 P-45.2.1 并列最优组。

`_p45_2_prefix_count` 数 `owned_atoms` 边界外的 claim 个数（局部 import L3 `iter_claims`，见 [[architecture/layer3-substituents]]）。

`_condensed_rank(info, cand)` 取**候选所辖缩合含氧酸中心的链内桥氧数**：`oxo_kind` 为 `phosphate` 或 `sulfate` 时，对 `covered_principal_ids` 内的 occurrence 取 `payload["oxo_z"]`（`_covered_oxo_z`），求 `_oxo_bridge_arms`（L1，中心 P/S 的 O-桥氧数）的最大值；**多核磷酸/硫酸以链中中心为功能母体（P-67.2.1）**，桥氧数多者优先；其余候选恒 0。

`_reorder_oxo_ester_side`（仅当候选全部 `oxo_kind == "phosphate"`）按 `_oxo_ester_side_score` 升序重排**缩合磷酸酯**的并列候选，把母体定向到糖/多元醇/甘油侧（P-67.1.3，数据驱动定向）。打分 `(臂上芳香含氮环数, -臂内氧数)` 越小越优先：`_oxo_ester_arm_carbons` 取中心 O-酯臂碳（中心-O-C，排除 `=O` 与桥氧），`_ester_arm_component` 沿「不越中心 P/S、不越与中心成键的酯氧/桥氧」求该臂的连通块，块内氧数取臂间最大、块所含芳杂环的吡啶型 N 单独计数。

出口契约：L2 只交 P-45.2.1 并列组；P-45.2.2 的位次集合须 L4 编号后才可得，由 L4 的 `prefix_locant_set` 给出并在 `namer._best_hit` 作末位裁决。`namer._candidate_phases` 取该组全部候选逐个跑 L3–L5。

### 关键数据类型

全部为冻结 dataclass 或 `str` 枚举，跨模块传递、可直接比较与哈希：

| 类型 | 定义位置 | 含义 |
|---|---|---|
| `SkeletonTopology` | `parent_skeleton` | 骨架拓扑：`ACYCLIC` / `RING_SYSTEM` |
| `ParentSkeleton` | `parent_skeleton` | 一个骨架候选：拓扑、原子序、覆盖的主基团 id |
| `SkeletonSelection` | `parent_skeleton` | 骨架候选集合（枚举结果或筛选后胜出者） |
| `PrincipalPriority` | `principal` | P-41 类号 + P-43 路径，`order=True` 可作比较键 |
| `PrincipalExpression` | `principal` | 主基团表达方式（当前只有 `SUFFIX`） |
| `PrincipalFeatureSpec` | `principal` | 主官能团规格：优先级、表达方式、锚点字段名 |
| `PrincipalGroupSelection` | `principal` | 选中的主官能团类及其全部 occurrence |
| `PrincipalParentSelection` | `parent_select` | 主官能团选择与骨架选择的冻结组合 |
| `PrincipalRelation` | `principal_expression` | 主基团与骨架的关系（骨架内/环外） |
| `PrincipalExpressionFacts` | `principal_expression` | 表达事实：类别、个数、关系、特征/锚点/附着原子 |
| `NumberingPolicy` / `ScaffoldSpec` / `ScaffoldIdentity` | `ring_scaffold` | 固定编号路径、骨架规格与身份 |
| `KindMeta` | `kind_registry` | kind 的中英词干、环型、环数、是否保留名 |
| `FusedNode` | `fused_system` | 稠环组分树节点（含稠合共享原子与附加组分） |
| `BridgeSegment` | `bridged_system` | 一段桥：两端桥头、桥内原子序与位次、上标位次对；`len()` 为桥内原子数 |
| `BridgedNode` | `bridged_system` | P-23 桥环母体：主环/主桥/次级桥划分 + 原子→位次映射 |
| `SpiroComponent` | `spiro_system` | 螺环系内一个单环组分（P-24.2）：环下标、原子集、螺原子 |
| `SpiroNode` | `spiro_system` | P-24.2 螺环母体：组分 + 螺描述符（含上标）+ 原子→位次映射 |
| `FbsNumbering` | `spiro_system` | 组分式螺环中一个组分的本位次编号候选（链序 + 位次标签） |
| `FbsComponent` | `spiro_system` | 一个环组分：`mono_ring`/`fused_ring`/`bridged_ring` + 基名/裸词干 + 编号候选 + 加氢前缀 |
| `FbsNode` | `spiro_system` | P-24.5~24.7 组分式螺环母体：引用序组分 + `links` + `cites` |

---

## 与其他层的契约

| 字段组 | 消费方 | 用途 |
|---|---|---|
| `owned_atoms` | L3 `claimable_block.iter_claims`、L2 自身 `_p45_2_prefix_count` | 取代基块的边界 |
| `chain` / `n_carbons` / `numbering_scaffold` / `scaffold_match` | L4 编号引擎 | 骨架原子序与固定编号映射 |
| `principal_expression_facts` | L3 位次来源、L5 后缀表达 | 附着原子/关系/电荷状态 |
| `stem_en`/`stem_zh` / `kind` | L5 词干引擎 | 母体名与后缀拼接 |
| `fused_tree` / `hydro_atoms` | L5 稠合组装、L4 加氢位次 | 稠环名与加氢前缀（`hydro_atoms` 亦由保留模板与生成式 HW 环给出） |
| `spiro_node` / `spiro_nodes` | L4 `spiro_numbering`、L5 `spiro_parent_names` | P-24.2 螺描述符候选裁决与 `spiro[a.b]` 词干 |
| `fbs_node` / `fbs_nodes` | L4 `fbs_numbering`、L5 `fbs_parent_name` | P-24.5~24.7 逐组分编号与组分式螺环名 |
| `bridged_node` / `bridged_nodes` | L4 `bridged_numbering`、L5 桥环词干 | P-23 编号候选裁决（P-14.4）与 von Baeyer 描述符/上标名 |
| `free_valence_order` | L5 `chain_engine` | 碳锚点自由价的 `-ylidene`/`-ylidyne` 词尾 |
| `radical_anchor_element` | L5 `assembler`、`assembler_prefixes` | 杂原子锚点自由基已并入母体名，故 L5 不另加前缀 |
| `single_atom_skeleton` | L5 | 单核母体阳离子/氮负离子/单核母体氢化物的收敛标记（骨架仅一个原子） |
| `anion` / `o_idx` / `thio_side` / `oxo_kind` / `n_oh` / `n_om` / `salt_meta` | L5 | 阴离子后缀、酯/硫代酯、含氧酸与盐组装 |
| `amide_z` / `hydrazide_n_idx` / `hydrazide_near_n_idx` | L5 `assembler`、`chain_engine`、`assembler_prefixes` | 硫代/亚氨酰胺词尾（P-43 类 16/17）与酰肼尾/肼基 N（P-66.3.1.1） |
| `heterane_z` / `lambda_n` / `lambda_atoms` / `stem_bare_en` / `stem_bare_zh` | L5 `assembler`、`chain_engine` | P-21 杂原子烃：氢化物词干、λ 标记与供链引擎拼 `ene`/`yne` 的裸词干 |
| `scaffold_id` / `scaffold_identity` | L4/L5 判读；`SPIRO_SCAFFOLDS` 被 L4 `indicated_hydrogen` 直接 import | 骨架身份标记（保留模板 id、`hw_mono`、`bridged`、`mono_spiro`/`fused_bridged_spiro`） |

`layer2` 的对外符号集中在 `parent_select`：`select_parent` 是唯一选母体入口，`finalize_parent_ownership` 供 `namer` 在候选补齐后固化所有权，`rule_driven_parent_candidates` 与 `select_principal_parent_skeletons` 可单独调用做骨架枚举。`ring_scaffold` 的模板表与 `_Q`/`_hydrogenated`、`standard_chain`、`extra_indicated_atoms`/`extra_hydrogenated_atoms`、`get_spec`、`component_stem`、`match_retained` 被 L4/L5 编号、指示氢与稠合组装模块直接引用；`spiro_system.SPIRO_SCAFFOLDS` 与 `hantzsch_widman` 的 `is_hw_scaffold`/`effective_valence`/`HW_RING_LAMBDA_Z` 是 L4 直接 import 的 L2 常量。扩展新环骨架的流程见 [[guides/adding-new-ring-system]]，跨层数据契约见 [[reference/core-data-contracts]]。

> **源:** `src/opensmiles/layer2/parent_select.py`、`principal.py`、`parent_skeleton.py`、`principal_expression.py`、`ring_scaffold.py`、`kind_registry.py`、`hantzsch_widman.py`、`fused_system.py`、`bridged_system.py`、`spiro_system.py`、`chain_walk.py`

相关页面：[[architecture/layer1-analyzer]]、[[architecture/layer3-substituents]]、[[architecture/layer4-numbering]]、[[architecture/layer5-name-assembly]]、[[concepts/atom-ownership]]、[[concepts/functional-group-priority]]、[[reference/core-data-contracts]]、[[guides/adding-new-ring-system]]。
