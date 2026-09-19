# Layer2: 母体选择 (Parent Selection)

> **源文件:** `src/namepredict/layer2/` 11 个 `.py`（含 `__init__.py`），约 2418 行 | **对外接口:** `parent_select.select_parent(info) -> list[dict]`、`parent_select.finalize_parent_ownership(parent, mol) -> dict`

---

## 概述

Layer2 接收 layer1 `analyze()` 的 info 字典，按 P-44 选出母体结构 (parent hydride)，返回 **P-45.2.1 并列最优候选组**；每个候选是一个 dict，供 layer3 提取取代基、layer4 在母体骨架上编号、layer5 组装名称。

| 文件 | 职责 |
|---|---|
| `parent_select.py` | 层唯一门面：P-44 编排、候选收集、原子归属、P-45.2 排序 |
| `principal.py` | P-41 注册表与主基团选择 |
| `parent_skeleton.py` | P-44 骨架枚举与筛选谓词 |
| `principal_expression.py` | 主基团表达 facts 与母体 dict 合成 |
| `ring_scaffold.py` | 骨架规格 + 保留 SMILES 模板表 + 环解析（唯一事实来源） |
| `kind_registry.py` | kind → `KindMeta`（中英词干/环型/环数/保留名） |
| `ring_expression_policy.py` | 环骨架上 typed 主基团表达的能力策略 |
| `fused_system.py` | P-25.3.2.4 稠环拆解为 `FusedNode` 树 |
| `bridged_system.py` | P-23 扩展 von Baeyer 桥环拆解为 `BridgedNode`（主环/主桥/次级桥 + 位次） |
| `chain_walk.py` | 脂肪族碳链行走（等长最长链全枚举） |

**输入 (info)**：`mol`（RDKit Mol）、`carbon_ids`/`n_carbons`、`fg_inventory`、`double_bonds`/`triple_bonds`、`rings`/`n_rings`/`has_ring`、`ring_systems`/`n_ring_systems`，以及 `namer._name_mol` 注入的 `root_ctx` 与 `salt`。FG 存在性一律由 `fg_inventory` 内容判定。

**输出 (parent dict)**：`kind`（母体类型）、`chain`（骨架原子序号）、`n_carbons`、`owned_atoms`（frozenset）、`covered_principal_ids`、`principal_occurrences`、`principal_group_count`、`principal_expression_facts`、`stem_en`/`stem_zh`、`numbering_scaffold`、`scaffold_id`/`scaffold_identity`/`scaffold_match`/`typed_ring_expression_supported`、`hydro_atoms`、`fused_tree` 或 `bridged_node`/`bridged_nodes`，以及按 FG 与拓扑分支写入的 `radical_c_idx`/`acyl_c_idx`、`double_bond(s)`/`triple_bond(s)`、`radical_ylidene`、`anion`、`o_idx`/`alkoxy_n`/`thio_side`、`hal_idx`/`hal_z`、`oxo_kind`/`n_oh`/`n_om`/`n_arms`/`salt_meta`。

**职责边界**：本层不切割取代基（L3）、不编号（L4）、不拼接名称字符串（L5），也不做烷烃兜底——无候选时返回空列表，由 `namer._run_candidates` 以 `no_assemblable_candidate` 显式失败。

层的调用链：

```
select_parent(info)
├─ _collect_candidates          rule_driven_parent_candidates → 剔除 None
│   ├─ select_principal_parent_skeletons   select_principal_group(P-41) + select_principal_skeletons(P-44)
│   └─ _express_selected        express_ring_principal / express_chain_principal
│       └─ _ring_scaffold_and_nodes  scaffold 身份 + P-25 稠环节点 + P-23 桥环节点（桥环接管身份）
├─ _finalize_ranked             pack_parent_stem(词干/编号 facts) → finalize_parent_ownership(owned_atoms)
└─ _reorder_p45_2(tied=True)    P-45.2.1 前缀取代基团数最大组
```

管线内部全部为无副作用纯函数（记忆化除外），可独立单测：L2 的输入只有 info，输出只有候选列表，不写回 info。

---

## 骨架候选（P-44）

`ParentSkeleton(topology, atom_ids, covered_principal_ids)` 与 `SkeletonSelection(candidates)`；`SkeletonTopology` 取 `ACYCLIC` / `RING_SYSTEM`。

`enumerate_principal_skeletons` = `_chain_candidates` + `_ring_candidates`：

- **链候选** `_open_chains`：以主基团锚点（`_anchors`）中的开环锚点为种子，`_all_chains_through` 全枚举等长最长链（平局交下游裁决），再加 `_pair_chains`（`_chain_through_two` 连接两锚点的链）；无锚点时退化为 `_longest_chain` 单候选（纯烃）。`_demoted_leaf_carbons` 把被压为前缀叶的中心碳列为禁行点（P-61.1.3）。同一原子集去重。
- **环候选** `_ring_candidates`：对 `info["ring_systems"]` 的每个环系构造候选；附着判据 `_ring_attaches`——锚点在环内即为附着，胺/醇/酮/自由基只认直接附着，其余类别允许经一个邻居间接附着。
- **覆盖度** `_chain_coverage`：链覆盖的 occurrence（胺任一臂在链即覆盖，P-62.2；其余要求全部锚点在链内）。

`select_principal_skeletons` 的筛选链：`keep_max_principal_coverage` → `keep_p44_1_2`（候选拓扑混合时用 `keep_senior_atom`，按 `_SENIOR_ATOMS` 判最优先元素）→ 全开链取 `keep_p44_3`（`p44_3_key`：杂原子数、原子数、元素计数），否则取 `keep_p44_2`（`p44_2_key`：有无杂原子、N 计数、senior 元素、完整环数、原子数、杂原子数）→ 末位 `keep_p44_4_unsaturation`（`p44_4_unsaturation_key` 计多重键数与双键数，经 `memo.by_key` 记忆；芳香键按半计）。

---

### P-44 谓词速查

| 谓词 / 键 | 条款 | 判据 |
|---|---|---|
| `keep_max_principal_coverage` | P-44 前置 | 保留覆盖主官能团最多的候选 |
| `keep_p44_1_2` / `keep_senior_atom` | P-44.1.1 / .1.2 | 候选拓扑混合时取含最优先元素者（`_SENIOR_ATOMS` 顺序 N>P>…>O>S>C） |
| `keep_p44_2` / `p44_2_key` | P-44.2 | 环候选：有杂原子 > N 计数 > senior 元素 > 完整环数 > 原子数 > 杂原子数 |
| `keep_p44_3` / `p44_3_key` | P-44.3 | 开链候选：杂原子数 > 原子数 > 元素计数元组 |
| `keep_p44_4_unsaturation` | P-44.4 | 多重键数 > 双键数（芳香键按半计） |

开链与环候选在枚举阶段并存；只要候选集中同时出现两种拓扑，就先过 P-44.1.2，随后按现存拓扑类型走对应的 P-44.2 或 P-44.3，最后统一过 P-44.4。

## 主基团选择（P-41）

`PRINCIPAL_REGISTRY` 由 `FG_SPECS` 中 `p41 != 0` 的条目派生，每项为 `PrincipalFeatureSpec(PrincipalPriority(p41_class, p43_path), expression, anchor_fields)`。

`select_principal_group(inventory)`：取非 `demoted` 条目对应的类别集合（`demoted` 叶 carboxy/cyano 不作主基团候选），过滤掉无规格者，按 `_effective_priority` 取 `min`；选中后取 `inventory.occurrences(选中类)`，得到 `PrincipalGroupSelection(group_class, occurrences)`；无候选时返回 `FG.NONE` 的空选择。

`_effective_priority(group_class, spec, inventory)` 给出候选类的实际 P-41 优先级：**含氧酸类（`FG.OXOACID`）的全部 occurrence 都满足 `oxoacid_is_acid(payload)`（酸式）时，升到 `OXO_ACID_P41`（类 7 内排在羧酸之后的等级）**，从而排在酯（9）与酰胺（11）之前；否则按注册表登记值参与比较。它是 `select_principal_group` 的比较键，也是「主基团互斥」的实现——只选单个最高优先级类别。

`_spec_from_fg` 只搬运注册表字段，`feature_spec` 为唯一查表入口；`PrincipalExpression` 只用到 `SUFFIX`。

注册表覆盖的 P-41 等级（`p41` 值）与主基团类的对应关系：自由基与酰基 1、羧酸 7、含氧酸 9（酸式时升为 `OXO_ACID_P41`）、酯 9、酰卤 10、磺酰胺与酰胺 11、腈 14、醛 15、酮 16、醇与硫醇 17、胺 19。`p43_path` 用于同级内的先后（如羧酸 `(1,)` 与含氧酸 `(0,)`、醇 `(1,)` 与硫醇 `(2,)`），`PrincipalPriority` 以 `order=True` 的冻结 dataclass 承载，可直接作 `min` 的比较键。

---

## 表达与 kind 正交化

两个入口 `express_chain_principal`（仅 `ACYCLIC`）与 `express_ring_principal`（仅 `RING_SYSTEM`），由 `parent_select._express_selected` 按骨架拓扑分派。

`express_ring_principal` 先调一次 `_ring_scaffold_and_nodes`，把 scaffold 身份与已算好的节点一次交给 `_ring_kind` 与 `_scaffold_fields` 复用。`_unsupported_typed_ring` 给出环酮被列为 typed 表达不支持的判据：主基团类为 `KETONE`、候选带 `scaffold_identity`、且 `typed_ring_expression_supported is False`。

`_parent_dict` 合成公共字段：`kind`、`chain`(=骨架原子序)、`n_carbons`、`covered_principal_ids`、`principal_occurrences`、`principal_group_count`、`principal_expression_facts`，与各类专属 fields 合并。

`PrincipalExpressionFacts`：`group_class`、`multiplicity`、`relation`（`IN_SKELETON`/`EXOCYCLIC`，按特征原子是否落在骨架内判定）、`occurrence_ids`、`characteristic_atoms`、`anchor_atoms`（官能团原锚点）、`attachment_atoms`（`_skeletal_attachments`：锚点全在骨架外时取其在骨架内的邻居）、`charge_state`（`NEUTRAL`/`ANION`/`MIXED`，由 `_is_anion_occurrence` 按周边原子形式电荷判定）。

### kind 决定

`_chain_kind(group_class, count, occurrences)`：

| 类别 | kind |
|---|---|
| `NONE` | `"none"`（仅 count == 0，纯烃） |
| `ACYL` | `"acyl"`（羰基头为 locant 1，P-65.1.7.2） |
| `RADICAL` | `"radical"` |
| `_OXO_FG_CLASSES`（`OXOACID`/`SULFONAMIDE`） | `_oxo_kind_of(occurrences)`：**取自 L1 payload 的 `oxo_kind`**，全部 occurrence 须同一 `oxo_kind`，否则返回 None 表示不支持 |
| 其余 | `group_class.value`（count ≥ 1） |

`_ring_kind` 复用同一张表并向 `_chain_kind` 透传 `occurrences`：自由基在无 scaffold 时返回 None（未知杂环无 -yl 词干，显式失败）；有 scaffold 且类别已注册时取 `_chain_kind`；环上外环 `-CHO` 单独放行 `"aldehyde"`（P-66.6.1.1.3）；否则 `_resolved_ring_kind` —— 苯与未注册稠环（`fused_hetero`）都收敛为 `"alkane"`，身份分别由 `scaffold_id` 与 `fused_tree` 承载；其余保留身份取 `scaffold.id` 本身（桥环即 `"bridged"`，身份由 `bridged_node` 承载）；`carbocycle` 或未命中模板时走 `_generic_ring_kind`，一律给 `"alkane"`（芳香多环与纯碳环的 kind 都收敛于此），身份仍由 `scaffold_id` 与 `fused_tree`/`bridged_node` 承载。

### 阴离子与酯/硫代酯字段

- `_ANION_FLAG_FGS = {ACID, OXOACID}`：该类别的全部 occurrence 都是阴离子时写 `anion: True`（`_expression_flags` 供链与环共用），L5 据此转 `-ate`。
- `_ester_o_idx(mol, payload)`：**羰基碳上另连烃基的单键 O 优先，无 O 时取 S**（硫代酯）。
- `_thio_fields(mol, o_idx)`：o_idx 为 S 时写 `{"thio_side": True}`（P-65.6.3.3.7.1，供 L5 换用 thioate 词尾与斜体 S 位次）。`ester_fields`（环骨架）与 `_chain_ester_fields`（链骨架）都补该标记；单 occurrence 时附 `o_idx` 与 `alkoxy_n: 0`。
- `_chain_acyl_halide_fields`：单一 occurrence 时从 `surr_idx` 取卤素，写 `hal_idx`/`hal_z`。

### 含氧酸字段与盐门控

`_chain_oxoacid_fields(info, occurrences, fields)`（单一 occurrence 才放行）写 `oxo_kind`/`n_oh`/`n_om`：

- kind **不在** `_OXO_Z_ANCHORED`（碳锚定的磺酸/膦酸等，中心不入母体）时直接返回字段、**不做盐门控**。
- kind 为中心自任母体者（磷酸/膦酸酯/硫酸酯）才校验金属盐：`n_om > 0` 时要求 `salt["n_metal"] == n_om`；`n_om == 0` 时不允许存在金属。通过后附 `n_arms` 与 `salt_meta`。
- 门控不通过返回 None，该候选被丢弃。

`express_chain_principal` 的含氧酸分派条件是 `selection.group_class in _OXO_FG_CLASSES`（含氧酸与磺酰胺共用）。

### 不饱和度与其他字段

`_chain_unsat_fields` 把骨架内 C=C/C≡C 写成 `double_bond`/`triple_bond`（单个）或 `double_bonds`/`triple_bonds`（多个）；`scaffold_id` 为 `carbocycle` 或 `bridged` 时由 `_kekule_ring_dbs` 补回被芳香感知剔除的环内 C=C，保留 mancude 母体覆盖的原子集（`_implied_ring_atoms`）则剔除其隐含多重键。自由基分支：`_mononuclear_radical` 把杂原子锚点收敛为表 2.1 单核氢化物骨架（S/N/P 的氧化态或键级并入词干），`_radical_ylidene` 判碳锚点自由价双键。`_semantic_anchor_fields` 为 `RADICAL`/`ACYL` 的单锚点写固定 locant 1 字段。

---

## kind 与词干注册

`kind_registry` 是母体 kind 的中英词干权威，注册表 `_REG` 只由 `ring_scaffold.all_specs()` 经 `_load_from_scaffold_specs()` 在 import 期填充：每个带 `stem_en`/`stem_zh` 的 `ScaffoldSpec` 生成一条 `KindMeta(kind, en, zh, ring, n_rings, retained)`。`get`/`parent_names` 是唯一读入口；未注册 kind 返回 None。

`pack_parent_stem(parent, mol)` 在 L2 出口补齐词干与编号字段：

- 按 `scaffold_id`（缺则 `kind`）查 `parent_names`；命中且尚未写 `stem_en`/`stem_zh` 时写入。
- 再由 `locant_prefix(scaffold_id)` 决定是否前置 locant 前缀（如 `1H-`、`1,3-`）；词干自带同一前缀时不重复前置；条件化前缀（`prefix_nh_conditional`）只在该环有未取代芳香 NH 时保留（`_ring_keeps_nh_prefix`），否则去前缀。
- 末步 `_attach_numbering_scaffold` 写入 `numbering_scaffold`（固定编号标签表），解析失败则原样返回。
- `mol` 缺省时把当前 mol 合并进候选，供词干判定取原子状态。

---

## 原子归属

`finalize_parent_ownership(parent, mol)` 生成不可变 `owned_atoms` frozenset（已是 frozenset 则原样返回），它是 L3 取代基块的边界：

- **`parent["kind"] in OXO_CENTER_KINDS`**（`phosphate`/`phosphonate`/`sulfate`，中心自任母体）：`owned_atoms` **只取主官能团特征原子、不含链**——链仅供编号，臂一律退为取代基。
- 其他 kind：`owned_atoms` = 链原子 ∪ 主官能团特征原子。

`_kind_fg_atoms` 的算法：

1. 种子 = occurrence 锚点 ∩ 链；种子为空（exocyclic，如苯甲酸的羧基）时改取与链相邻的锚点。
2. **`facts.group_class in _WHOLE_FG_ATOMS = {OXOACID, SULFONAMIDE}`** 时，只并「本骨架覆盖的 occurrence」（`covered_principal_ids`）的 `characteristic_atoms` 后返回——**不沿线借臂**，因为中心 P/S 不是碳骨架成员、没有可外借的臂。
3. 其他类别从种子出发，沿邻接并入属于 occurrence 特征原子、且本身不是锚点的邻居（沿线借臂）。

---

## P-45.2 排序

`select_parent(info)` = `_collect_candidates` → `_finalize_ranked` → `_reorder_p45_2(tied=True)`：

- `_collect_candidates`：`rule_driven_parent_candidates`（= `select_principal_parent_skeletons` + `_express_selected`）剔除 None。
- `_finalize_ranked`：逐个候选先 `pack_parent_stem(c, mol)`（补 `stem_en`/`stem_zh`、locant 前缀与 `numbering_scaffold`），再 `finalize_parent_ownership` 固化 `owned_atoms`。
- `_reorder_p45_2`：候选不足 2 个直接返回；排序键为三元组 **`(前缀取代基团数, 缩合磷酸链内桥氧数, 原始序)`**，前两项降序、原始序升序稳定排序；`tied=True`（层出口）只保留第一项并列最大的一组，即 P-45.2.1 并列最优组。

`_p45_2_prefix_count` 数 `owned_atoms` 边界外的 claim 个数（局部 import L3 `iter_claims`）。

`_condensed_rank(info, cand)` 取**候选所辖缩合含氧酸中心的链内桥氧数**：`oxo_kind` 为 `phosphate` 或 `sulfate` 时，对 `covered_principal_ids` 内的 occurrence 取 `payload["oxo_z"]`，求 `_oxo_bridge_arms`（L1，中心 P/S 的 O-桥氧数）的最大值；**多核磷酸/硫酸以链中中心为功能母体（P-67.2.1）**，桥氧数多者优先；其余候选恒 0。

出口契约：L2 只交 P-45.2.1 并列组；P-45.2.2 的位次集合须 L4 编号后才可得，由 L4 的 `prefix_locant_set` 给出并在 `namer._best_hit` 作末位裁决。`namer._candidate_phases` 取该组全部候选逐个跑 L3–L5。

---

## 环骨架识别

`ring_scaffold.py` 是骨架身份的唯一事实来源。

**模板表**：`_TEMPLATES` 共 **88 条**保留 SMILES 模板，每条含 `smiles`/`naming_class`/`stem_en`/`stem_zh`/`fused`，另有可选 `fused_prefix`、`fused_stem`、`locant_prefix`、`prefix_nh_conditional`、`standard`。命名类分布：`monohetero` 55、`naph_family` 9、`fused56` 8、`xanthene` 2，以及 `mono_carbo`/`anthra`/`phenanthrene`/`pyrene`/`picene`/`chrysene`/`carbazole`/`acridine`/`phenothiazine`/`benzodioxole`/`purine`/`indolizine`/`pyrrolizine`/`steroid` 各 1。饱和五元双杂环（如 `pyrazolidine`、`thiadiazolidine124`）取 Hantzsch-Widman 名而非氢化芳环名。

`_spec_from_template` 由模板派生 `ScaffoldSpec`（环数与环型由查询分子自动算、`retained=True`、`numbering.standard_path` 取 `standard` 标签），`all_specs()`/`get_spec()` 对外供全库使用；`kind_registry._load_from_scaffold_specs()` 在 import 期据此建 `KindMeta` 词干权威，`pack_parent_stem` 按 `scaffold_id`（缺则 `kind`）取母体名。

**固定编号**：41 条模板携带 `standard = (labels, order)`，含 `dihydrofuran`/`dihydropyran`/`dihydropyrrole`/`dihydroimidazole`/`dihydrothiazole` 这类部分不饱和环（字面位次即固定编号，P-14.4(a)/(b)）；import 期 `_validate_standard_fields` 校验 `order` 是 `0..n-1` 的排列、`labels` 长度等于模板原子数。派生 `_STANDARD_LABELS`/`_STANDARD_ORDERS`；`standard_chain` 把固定编号映射为分子原子序，`numbering_scaffold_facts` 在标签数与骨架原子数一致时物化 `numbering_scaffold`。

**匹配与身份**：`match_retained`/`_match_with_map` 先按元素签名（`_elem_sig`/`_TEMPLATE_ELEM`）剪枝，再要求子图同构原子集**精确等于**骨架原子集，并由 `_is_induced_match` 校验原子集诱导子图的键数等于模板键数——子图同构容忍目标多出的键，缺此校验则环系的真子图模板会把多出的环静默丢掉。`_isolated_saturated_ring` 判原子集是否为不与他环稠合、且 Kekulé 视图环内全为单键的单环（RDKit 可能误判其芳香），是则只允许匹配全单键模板（P-31.2）；`mancude_only` 时跳过非 `fused` 的饱和保留名（不作稠合零件）。精确匹配失败后按 `_hydrogenated` 全氢化骨架再比对（`memo.by_mol` 记忆，P-25.3.4），同样施加诱导覆盖与全单键两重校验。`resolve_ring_scaffold` 命中则返回 `spec.identity`，否则 `_generic_carbocycle`：按骨架内完整环数返回 `fused_hetero`（≥ 2）/ `carbocycle`（单环，`ring` 记 `carbo`）。

**命名类速查**：

| naming_class | 条数 | 代表模板 | 环表达策略 |
|---|---|---|---|
| `monohetero` | 55 | 呋喃/吡啶/吡喃/噁唑/哌啶/氧杂环己烷 | 仅环内酮（多类）、无醇/胺白名单 |
| `naph_family` | 9 | 萘、喹啉、异喹啉、色烯、蝶啶 | 环内醇/酮 |
| `fused56` | 8 | 吲哚、苯并呋喃、苯并噻唑、苯并二氧杂环戊烯 | 环内酮 |
| `xanthene` | 2 | 呫吨、噻吨 | 环内酮、环内醇 |
| `mono_carbo` | 1 | 苯（kind 收敛 `alkane`） | 环内醇/酮/胺、环外酸 |
| `carbocycle`（非模板） | — | 未注册单环（无模板命中时按环数兜底） | 环内醇/酮/胺、环外酸 |
| `fused_hetero` / `fused`（非模板） | — | 未注册稠环 | 环内醇/酮/胺、环外酸/腈 |
| `bridged`（非模板） | — | 未注册桥环 | 环内醇/酮/胺、环外酸/腈 |
| 其他单例 | 各 1 | 蒽/菲/芘/苉/屈/咔唑/吖啶/吩噻嗪/嘌呤/中氮茚/吡咯嗪/甾体 | 环内酮（`purine`/`carbazole`/`acridine`/`phenothiazine` 并入单杂环组）；`picene`/`indolizine`/`pyrrolizine` 未列入白名单 |

**加氢与指示氢**：`hydrogenated_atoms`（模板 Kekulé 双键位在分子中已饱和者；环杂原子上的 H 交由指示氢承载；计数落在 `HYDRO_MULT_N` 之外且存在环外多重键时，其中一位由指示氢承载，P-31.2.2/P-58.2.1）、`extra_indicated_atoms`（P-58.2.1）、`mancude_ring_atoms`（保留 mancude 母体覆盖的整个不饱和环）、`_kekule_double_atoms`（模板双键端点缓存）。`locant_prefix` 返回 `(en, zh, nh_conditional)`，条件化前缀（如 `1H-`）只在环内存在未取代芳香 NH 时保留（`_ring_keeps_nh_prefix`）。

**稠合前缀**：`retained_fusion_prefix` 先查模板 `fused_prefix`，否则回退 `fusion_carbocycle_prefix`；`component_stem` 给出组分的稠合词干（`fused_stem` 优先，否则经 `_bracketed_stem`）。要点：

- 饱和环作稠合组分须取 **mancude 对应名**（P-25.3.1.2.1 / .2.3），而不是把饱和环名去尾加「并」：`oxane → ("pyrano", "吡喃并")`、`thiolane → ("thieno", "噻吩并")`、`thiane → ("thiopyrano", "噻喃并")`；`azepane` 的 `fused_stem` 为 `("azepine", "氮杂卓")`（P-25.3.1.2.2）。
- `_bracketed_stem`：杂原子位次前缀为纯数字（`1-`/`1,3-`/`1,2,4-`）且词干本身不以数字开头时，稠合组分名取方括号形式（`[1,3]oxazole`、`[1]benzofuran`、`[1,3]benzodioxole`），即 P-25.3.5 的「稠合母体须引用完整位次」；`1H-` 型前缀与词干已带位次者（`1,2-oxazole`）不加方括号。
- `_FUSION_CARBOCYCLES`（6 条单环烃附加组分，P-25.3.2.2.1）：`cyclohexane` 的稠合前缀是保留前缀 **`benzo`/`苯并`**（该条「除 benzo 外」的例外），其余为 `cyclopropa`/`cyclobuta`/`cyclopenta`/`cyclohepta`/`cycloocta`。`match_fusion_carbocycle` 做精确匹配，`omits_fusion_numbers` 判 benzene 与单环烃附加组分省略数字位次（P-25.3.8.1）。

---

## 稠环拆解与多环母体

`_ring_scaffold_and_nodes` 取与骨架原子集一致的环系 dict（`_skeleton_system`），其 `sssr_indices` ≥ 2 时先调 `decompose_fused_system(info, system)` 得 P-25 稠环节点，再调 `try_bridged_scaffold` 得 P-23 桥环节点；桥环命中时以 `bridged[0].scaffold_identity()` 接管 scaffold 身份，`kind` 由 `_resolved_ring_kind` 给出 `"bridged"`。两套拆解都独立于 scaffold 模板：未注册系统同样产出。

`_scaffold_fields(info, skeleton, facts, scaffold, fused_tree, bridged)` 按传入节点写入字段：`bridged` 非空时写 `bridged_node`（并列最优的首个候选）与 `bridged_nodes`（全部并列候选，供 L4 裁决），否则 `fused_tree` 非空时写 `fused_tree`——两者互斥，桥环优先。`scaffold=None` 时该函数自行调 `_ring_scaffold_and_nodes`。

`FusedNode`：`scaffold_id`、`atom_ids`、`ring_indices`、`fusion_shared`（与母体的稠合共享原子）、`attached`（附加组分，递归）、`fused_stem`、`fused_prefix`、`fused_omit_numbers`。

拆解流程：`_candidates_for` 从每个自身可作稠合零件的环种子出发，沿融合图（`_fusion_adj`）增长枚举候选环集，用元素计数超集剪枝（`_has_template_superset`）与原子集去重；`_select_base` 按 P-25.3.2.4 准则 (a)–(j) 逐条 `narrow` 收窄——(a) 最优先杂原子、(b) 环数多、(c) 环大小降序、(d) 杂原子总数、(e) 杂原子种类、(f) 最优先杂原子计数、(g) 水平行环数、(h) 杂原子位次最低，逐元素位次最低，稠合碳位次最低；准则 (g)–(j) 依赖 L4（`preferred_orientations` + `number_fused_system`，经 `_numbered_locants`）算位次。全部并列时取环集升序最小。`_decompose` 递归：选出母体组分后，剩余环按融合图连通分量（`_ring_components`）递归为附加组分，共享原子记入 `fusion_shared`。

**环表达策略**：`ring_expression_policy._POLICIES` 是 (命名类, 主基团类, 关系) 三元组白名单，`supports_ring_expression` 判定 `typed_ring_expression_supported`。白名单覆盖 carbocycle/mono_carbo/naph_family 的环内醇/酮/胺与环外酸，单杂环与稠杂环（`monohetero`/`fused56`/`purine`/`carbazole`/`acridine`/`phenothiazine`/`benzodioxole`/`anthra`/`phenanthrene`/`pyrene`）、`xanthene`/`steroid`、未注册稠环（`fused_hetero`/`fused`）与桥环（`bridged`）的环内酮/醇/胺与环外酸/腈。该判据与 `scaffold_id`/`scaffold_identity`/`scaffold_match`/`hydro_atoms`、`fused_tree`/`bridged_node` 一并写入 parent，供 L4 取固定编号与加氢位次。

## 桥环拆解（P-23）

`bridged_system.py` 表达 P-23 扩展 von Baeyer 桥环：把环系拆成主环、主桥与次级桥，并给出全原子位次。

**数据结构**：`BridgeSegment(heads, atoms, numbers, locants)` 描述一段桥，`heads` 为两端桥头、`atoms` 沿 `heads[0] → heads[1]` 序；`__len__` 返回桥内原子数，即 von Baeyer 描述符中的数字（P-23.2.6.1.2），`reversed()` 掉转方向。`BridgedNode(atom_ids, numbering, main_ring, ring_segments, main_bridge, independent_bridges, dependent_bridges, ring, n_rings)` 是桥环母体：`main_ring` 为编号序的主环原子、`ring_segments` 为两条环段的长度、`ring` 取 `carbo`/`hetero`。只读派生属性：`secondary_bridges`（独立桥在前、依赖桥在后，即引用顺序）、`descriptor`（主环两段长 + 主桥长 + 各次级桥长）、`locant_pairs`（次级桥的上标位次对，小者在前，按引用顺序）；`scaffold_identity()` 返回 `ScaffoldIdentity("bridged", "bridged", n_rings, ring)`，id 与命名类同名。

**识别流程**：`decompose_bridged_system(info, system)` 是公共入口，非桥环返回空表。前置筛：原子数 ≥ 4、`_ring_count`（键数 − 原子数 + 1，P-23.2.6.1.1）算出的环数 ≥ 2、桥头 `heads` 取度 ≥ 3 的原子（P-23.1.1）且至少两个。`_bridges` 拆出全部桥：每个「去桥头分量」须为简单路径且两端各邻接一个不同桥头（`_path_segment`，单原子分量被两桥头夹持时即为 1 原子桥），再加桥头之间的直键（0 原子桥，P-23.1.2）；并要求 `len(bridges) - len(heads) + 1` 等于环数。`_candidates` 在桥图 H（顶点 = 桥头、边 = 桥，`_h_graph`）上枚举简单环（`_cycles`，上限 `_MAX_CYCLES`，两顶点间的平行边视为二元环），以环上边集为主环、其中一段为主桥、该段两端为主桥头对；其余桥须构成以主环桥头为端点的简单链（`_secondary_ok`：主环外顶点度恰为 2 且无环）。`_to_node` 给候选编号，编号失败返回 None；候选组装后还要求描述符数字之和加 2 等于骨架原子数（P-23.2.6.1.4），不符者丢弃。

**收窄阶梯** `narrow_candidates` 依次施加（每步走 `layer4.numbering_engine.narrow`，末步不做确定性 tie-break）：主环原子数最大（P-23.2.1）→ 主桥最长（P-23.2.4）→ 主环两段长之差最小（P-23.2.6.2.1）→ 独立桥长度序（P-23.2.6.2.2）→ 依赖桥最少（P-23.2.6.2.3）→ 上标位次集合最低（P-23.2.6.2.4）→ 引用序位次最低（P-23.2.6.2.5）。

**编号**：主环从长段一端进入——长边在前（P-23.2.3），必要时把两条环段的行进方向一并掉转以保持路径连续；主环依次占位次 1..n，主桥自桥头续编，`_number_one` 从紧邻较高编号桥头的一端开始给桥内原子编号（P-23.2.6.3），`locants` 记两端桥头的位次对。次级桥按引用顺序编号：独立桥（两端桥头都已入编号，即在主环或主桥上）在前，依赖桥随后，`_number_secondary` 逐个取当前可编号者，出现不可编号者则该候选失败。引用顺序由 P-23.2.6.2.2/.5 裁决：`_orders` 固定长桥在前、只对同长度桥排列，枚举数超过 `_MAX_ORDERS` 时退化为「长度降序」单一顺序；逐个顺序编号后取（次级桥长度降序、引用序位次展平）键最小者。

**路由** `try_bridged_scaffold(info, scaffold, fused_tree, system)`：先取 `decompose_bridged_system`，空表即让位。**条件 A** —— `_fusion_naming_applies` 为假时直接走桥环，该判据为「环系内至少两个 ≥ 5 元环」（P-52.2.4.1）且 `_cannot_be_mancude` 为假（环系内已成满 4 根 σ 键的碳腾不出 π 键，写不出最大非累积双键，P-25.3.1.2）。否则要求 scaffold 身份为 `fused_hetero`（身份为 None 或是保留模板时维持既有稠环路径、不走桥环）、`_tree_covers_rings` 确认稠环树的 `ring_indices` 覆盖环系的全部 `sssr_indices`（盖不住环系的全部环则回退桥环）；末步 `_p25_names_ok` 反问 `layer5.fused_namer.fused_parent_names`——**能拼出稠合词干（条件 B）则保持稠环路径**，拼不出才回退桥环。

## 关键数据类型

全部为冻结 dataclass 或 `str` 枚举，跨模块传递、可直接比较与哈希：

| 类型 | 定义位置 | 含义 |
|---|---|---|
| `SkeletonTopology` | `parent_skeleton` | 骨架拓扑：`ACYCLIC` / `RING_SYSTEM` |
| `ParentSkeleton` | `parent_skeleton` | 一个骨架候选：拓扑、原子序、覆盖的主基团 id |
| `SkeletonSelection` | `parent_skeleton` | 骨架候选集合（枚举结果或筛选后胜出者） |
| `PrincipalPriority` | `principal` | P-41 类号 + P-43 路径，`order=True` 可作比较键 |
| `PrincipalFeatureSpec` | `principal` | 主官能团规格：优先级、表达方式、锚点字段名 |
| `PrincipalGroupSelection` | `principal` | 选中的主官能团类及其全部 occurrence |
| `PrincipalParentSelection` | `parent_select` | 主官能团选择与骨架选择的冻结组合 |
| `PrincipalRelation` / `PrincipalChargeState` | `principal_expression` | 主基团与骨架的关系、电荷状态 |
| `PrincipalExpressionFacts` | `principal_expression` | 表达事实：类别、个数、关系、特征/锚点/附着原子、电荷态 |
| `NumberingPolicy` / `ScaffoldSpec` / `ScaffoldIdentity` | `ring_scaffold` | 固定编号路径、骨架规格与身份 |
| `KindMeta` | `kind_registry` | kind 的中英词干、环型、环数、是否保留名 |
| `RingExpressionPolicy` | `ring_expression_policy` | 一条能力策略：命名类集合 + 主基团类 + 允许的关系 |
| `FusedNode` | `fused_system` | 稠环组分树节点（含稠合共享原子与附加组分） |
| `BridgeSegment` | `bridged_system` | 一段桥：两端桥头、桥内原子序与位次、上标位次对；`len()` 为桥内原子数 |
| `BridgedNode` | `bridged_system` | P-23 桥环母体：主环/主桥/次级桥划分 + 原子→位次映射 |

---

## 与上下游层的契约

| 字段组 | 消费方 | 用途 |
|---|---|---|
| `owned_atoms` | L3 `claimable_block.iter_claims`、L2 自身 `_p45_2_prefix_count` | 取代基块的边界 |
| `chain` / `n_carbons` / `numbering_scaffold` / `scaffold_match` | L4 编号引擎 | 骨架原子序与固定编号映射 |
| `principal_expression_facts` | L3 位次来源、L5 后缀表达 | 附着原子/关系/电荷态 |
| `stem_en`/`stem_zh` / `kind` | L5 词干引擎 | 母体名与后缀拼接 |
| `fused_tree` / `hydro_atoms` | L5 稠合组装、L4 加氢位次 | 稠环名与加氢前缀 |
| `bridged_node` / `bridged_nodes` | L4 `bridged_numbering`、L5 桥环词干 | P-23 编号候选裁决（P-14.4）与 von Baeyer 描述符/上标名 |
| `anion` / `o_idx` / `thio_side` / `oxo_kind` / `n_oh` / `n_om` / `salt_meta` | L5 | 阴离子后缀、酯/硫代酯、含氧酸与盐组装 |

`layer2` 的对外符号集中在 `parent_select`：`select_parent` 是唯一选母体入口，`finalize_parent_ownership` 供 `namer` 在候选补齐后固化所有权，`rule_driven_parent_candidates` 与 `select_principal_parent_skeletons` 可单独调用做骨架枚举。`ring_scaffold` 的模板表、`_Q`/`_hydrogenated`、`standard_chain`、`extra_indicated_atoms`、`get_spec` 被 L4 编号与指示氢模块直接引用。

> **源:** `src/namepredict/layer2/parent_select.py`、`principal.py`、`parent_skeleton.py`、`principal_expression.py`、`ring_scaffold.py`、`kind_registry.py`、`ring_expression_policy.py`、`fused_system.py`、`bridged_system.py`、`chain_walk.py`
