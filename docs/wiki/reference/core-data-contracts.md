# 核心数据合约 (Core Data Contracts)

> NamePredict 六层命名流水线跨层数据结构的字段规格；字段名与语义以 `src/namepredict` 源码为唯一依据。

## 端到端数据流图

```text
SMILES ─preprocess─► Mol（同位素采集 isoNuclide/iso2H/iso3H · 互变异构/酸电荷归一）
  │ _simple_species 查表命中 → 直接给保留名（SIMPLE_MOLECULES）
  │ dissociate_salt → (organic Mol, salt)
  │ analyze(organic) → info dict（12 键，含 fg_inventory / lambda_atoms）
  │ namer._name_mol 注入 root_ctx、salt
  │ select_parent(info)  （preset → _finalize_ranked → P-45.2 / P-67 重排）
  ▼
list[parent dict] ─finalize_parent_ownership─► owned_atoms
  │ extract_substituents(info, parent)
  │   ├ iter_claims → ClaimedBlock      (layer3/claimable_block)
  │   └ SubstituentNamer.name → SubstituentName
  ▼
list[subst dict] ─number(parent, subs)─► numbered dict
  │   └ 桥环：bridged_nodes(L2 候选) → bridged_node(L4 写回) → L5 描述符/上标
  │ assemble(numbered)  （取名 → hydro/环鎓/核素 → 环 oxide → 前缀 → 阴离子/R-S/盐）
  ▼
NameResult ─namer._ok_result─► name(smiles)
```

组装失败或分子多片段时，`namer._name_mol` 回落到 `_name_components`：唯一可配对的阴阳离子走 `pair_unique_ions` 的 P-77 二元盐名，其余组分名字母序以 `COMPONENT_SEP = "; "` 连接，重复组分子 `_merge_repeats` 加 `bis`/`tris`。

覆盖门控 `build_coverage_ledger(mol, *, owned_atoms, names)` 只判 gap/overlap（`namer._prepare_candidate` 链路）。

---

## 1. info dict（L1 输出）

产出 `analyzer.analyze(mol)`（`base` + `_collect_fgs` + `_ring_meta` + `lambda_atoms`，共 12 键），`namer._name_mol` 另注入 `root_ctx`、`salt`；清单取用入口 `inventory_from_info`（缺失抛 `KeyError`）。

| 字段 | 类型 | 含义 | 生产 → 消费 |
|---|---|---|---|
| `mol` | `Mol` | 去盐后的有机部分 | `analyze` → 各层 |
| `carbon_ids` / `n_carbons` | `list[int]` / `int` | 碳原子索引及数量 | `analyze` → L1/L2 骨架枚举 |
| `double_bonds` / `triple_bonds` | `list[dict]` | `{"c1","c2"}`，C=C / C≡C（非芳香） | `_cc_bond_entries` → L2 `_chain_polys` |
| `fg_inventory` | `FunctionalGroupInventory` | 官能团事实的唯一出口 | `build_inventory` → L2/L3 |
| `rings` / `n_rings` / `has_ring` | `list[dict]` / `int` / `bool` | SSSR 环及数量、是否含环 | `_ring_meta` → `namer._prepare_candidate` |
| `ring_systems` / `n_ring_systems` | `list[dict]` / `int` | 稠合 + 螺环合并的环系（`atom_ids`/`sssr_indices`/`fusion_edges` 等） | `_ring_meta` → L2 `_scaffold_fields`、L4 |
| `lambda_atoms` | `dict[int, int]` | `{原子 idx: 键数}`，只收电中性且键数偏离标准值者（P-14.1.2） | `tools.lambda_notation.nonstandard_bonding` → L2 杂原子链 |
| `root_ctx` | `tuple[Mol, list[int]]` | `(根分子, 本分子原子→根索引)`；顶层为 `(organic, range(n))` | `namer._name_mol` → L3 在根分子上重算 R/S |
| `salt` | `dict` | L0 盐元数据（见下） | `namer._name_mol` → L2 盐门控（写 `salt_meta`）、L5 盐后缀 |

### L0 盐元数据（`salt` dict）

`dissociate_salt` 调 `_from_frags`：逐片段经 `_frag_role(mol) -> (role, en, zh) | None` 定角色（`metal` 阳离子 / `halide` 阴离子 / `hx` 中性卤化氢，None 即有机片段），按存在的角色分支产元数据；不满足简单盐条件时返回 `(mol, {})`。

| 键 | 类型 | 含义 |
|---|---|---|
| `metal` / `metal_zh` / `n_metal` | `str` / `str` / `int` | 金属阳离子（`METAL_ION_EN`/`METAL_ION_ZH`）；多原子阳离子 `[NH4+]`（`azanium`/铵）经 `POLY_CATION` 并入此通道 |
| `halide` / `halide_zh` | `str` / `str` | 卤素阴离子 X⁻（`HALIDE_EN`/`HALIDE_ZH`）；多原子阴离子（`hydroxide`/`nitrate`/`perchlorate`/`azide`/`hexafluorophosphate`）经 `POLY_ANION` 并入此通道，多份用 `MULT_*` 加数量前缀 |
| `acid_salt` / `acid_salt_zh` | `str` / `str` | 氢卤酸加合物 HX（`HALIDE_HX_*`，P-71.3）；中性有机碱 + 单一 Cl⁻ 亦归此键 |
| `n_org` | `int` | 有机片段份数（≥1）；金属盐臂 >1 时 L5 用 `bis(...)` 括起 |

非盐专用件：`pair_unique_ions(frags) -> (阳离子片段, 阴离子片段, n_阳, n_阴) | None`（唯一配对且整体电荷守恒，供 P-77 二元盐名）；`constants.SIMPLE_MOLECULES` / `SIMPLE_MAX_HEAVY` / `ELEMENT_METAL_NAMES` 供 `namer._simple_species` 做重原子 ≤ 3 的单质/简单分子查表（命中直接产名，`meta.reason = "simple_species"`）。

---

## 2. FG occurrence / inventory

`FunctionalGroupOccurrence`（`layer1/functional_group_inventory.py`）：

| 字段 | 类型 | 含义 |
|---|---|---|
| `id` | `str` | `"{fg}:{i}"`，如 `"acid:0"` |
| `group_class` | `FunctionalGroupClass` | 18 个注册 FG 类别（含 `AZANIDE`/`THIONE`/`HETERANE`）+ `NONE = "alkane"` |
| `characteristic_atoms` | `frozenset[int]` | 默认 `center_idx` ∪ `surr_idx`；例外见 `FG_ATOM_FNS` |
| `parent_anchors` | `frozenset[int]` | 由 `FgSpec.anchors` 声明的 key 从 payload 取出 |
| `payload` | `dict` | L1 原始条目 |
| `demoted` | `bool` | P-41 降级前缀叶（carboxy/cyano），碳排除出主链 |

`FunctionalGroupInventory`：字段 `entries`、`has_anion`（分子内是否存在负形式电荷原子）；方法 `occurrences(group_class)`（该类未降级条目）、`demoted_entries()`；构建入口 `build_inventory(lists, mol, demoted, has_anion)`。

**P-41 仲裁**（`_arbitrate_parts`，先于清单构建）：先 `_drop_mixed_anion_os`/`_drop_mixed_anion_acids` 让同族阴离子条目胜出；`_SUPPRESSIBLE`（acid/ester/acyl_halide/amide/nitrile/aldehyde）遇 `p41` 更小者整组退出；`_LEAF_DEMOTED`（acid/nitrile）标 `demoted`、碳排除出主链，其余组合 FG 整组清空；`_PRESENCE_SKIP`（oxoacid/sulfonamide）不判存在性；`has_anion` 时把 `cation` 从 present 中剔除（阴离子 > 阳离子）。返回 `(parts, demoted_ids)`。

**payload 键**：

| 类别 | payload 键 |
|---|---|
| `radical` | `center_idx`、`surr_idx = []` |
| `acyl`、`acid`、`ester`、`acyl_halide`、`amide`、`nitrile`、`aldehyde`、`ketone`、`thione`、`alcohol`、`thiol`、`amine`、`azanide`、`cation`、`heterane` | `center_idx`、`surr_idx` |
| `oxoacid`、`sulfonamide` | `oxo_z`、`center_idx`、`oxo_kind`、`n_oh`、`n_om`；**仅碳锚定条目**另有 `surr_idx`，只收阴离子氧（中心自任母体时无该键） |

- `anion_os`：仅 `alcohol`/`thiol`，中心 O⁻/S⁻ 为真（`_mark_anion_os`，P-41 类 4 定级）。
- `n_oh` = 中性含 H 氧数；`n_om` = 阴离子氧数；`n_arms` 不在 payload，L2 取 `payload.get("n_arms", 0)`。
- `oxo_kind`：B — `boronic`；P — `phosphate`/`phosphonate`；S — `sulfate`/`sulfonic`/`sulfonate`/`sulfonamide`/`sulfonyl_chloride`；`_OXO_CLASS_BY_KIND` 只分走 `sulfonamide`。
- 锚点 key（`FG_SPECS.anchors`）：`center_idx` 用于 radical/acyl/azanide/cation/acid/oxoacid/sulfonamide/ester/acyl_halide/amide/nitrile/aldehyde/ketone/thione/heterane；`alcohol`/`thiol`/`amine` 取 `surr_idx`。
- 特征原子例外 `FG_ATOM_FNS`：`oxoacid`/`sulfonamide` → `_oxoacid_atoms`（中心 P/S + 非碳邻居 + 锚点碳）；`azanide`/`heterane`/`cation` → `_cation_atoms`（只中心）；`nitrile` → `_nitrile_atoms`（只腈碳与三键氮）。

SMARTS 表 `layer1/fg_local_smarts.FG_SMARTS` 共 202 条、17 个 FG 键（heterane 一族按元素批量生成）。

---

## 3. parent dict（L2 输出）

产出 `select_parent(info)` → `list[dict]`：`_collect_candidates`（`rule_driven_parent_candidates` → `express_ring_principal`/`express_chain_principal` → `_parent_dict`）→ `_finalize_ranked`（`pack_parent_stem` + `finalize_parent_ownership`）→ `_reorder_p45_2(..., tied=True)` → `_reorder_oxo_ester_side`。整组并列最优逐个尝试 L4–L5。

| 字段 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `kind` | `str` | FG 类别名 / `oxo_kind` / `"radical"` / `"acyl"` / `"cation"` / `"azanide"` / `"heterane"` / 纯烃 `"alkane"` / 桥环 `"bridged"` / 生成式 `"hw_mono"` | `_chain_kind` / `_ring_kind` → L4/L5 |
| `chain` / `n_carbons` | `list[int]` / `int` | 骨架原子顺序及长度 | `_parent_dict`；L4 `number()` 覆盖为编号顺序 |
| `covered_principal_ids` | `tuple[str,...]` | 本骨架表达的主基团 id | `_parent_dict` → L2 |
| `principal_occurrences` | `tuple[FunctionalGroupOccurrence,...]` | 主基团 occurrence 全集 | `_parent_dict` → L2 |
| `principal_group_count` | `int` | 本骨架覆盖数 | `_parent_dict` → L2 |
| `principal_expression_facts` | `PrincipalExpressionFacts` | 主基团表达事实（见下） | `_facts` → L3/L4 |
| `radical_c_idx` / `acyl_c_idx` | `int` | 固定 locant 1 单锚点字段（仅 RADICAL/ACYL） | `_semantic_anchor_fields` → L4 |
| `anion` | `bool` | 酸全阴离子（仅 ACID/OXOACID/ALCOHOL/THIOL） | `_expression_flags` → L5 |
| `oxo_kind` | `str` | 含氧酸 kind | `_chain_oxoacid_fields` → L5 |
| `n_oh` / `n_om` / `n_arms` | `int` | 酸式 H 氧 / 阴离子氧 / O–R 臂数 | 同上 → L5 |
| `salt_meta` | `dict`/`None` | 盐门控后的 L0 盐元数据；`n_om > 0` 须 `n_metal == n_om` | 同上 → L5 |
| `o_idx` / `alkoxy_n` | `int` | 酯侧杂原子索引（硫酯为 S）/ 臂计数；`kind in ("ester","acid")` 均写 | `_chain_ester_fields` / `ester_fields` → L3/L5 |
| `thio_side` | `bool` | 硫代酯标记（酯氧换 S） | `_thio_fields` → L5 |
| `hal_idx` / `hal_z` | `int` | 酰卤卤原子索引与原子序 | `_chain_acyl_halide_fields` → L5 |
| `amide_z` | `int` | 酰胺双键杂原子 Z（`N` 亚氨酰胺 / `S` 硫代酰胺，P-43 类 16/17） | `_amide_fields` → L5 |
| `hydrazide_n_idx` / `hydrazide_near_n_idx` | `int` | 酰肼远 N（N′）/ 近 N（P-66.3.0） | `_amide_fields` → L3/L5 |
| `radical_anchor_element` | `str` | 杂原子锚点自由基的元素符号 | `_mononuclear_radical` → L5 |
| `free_valence_order` | `int` | 碳锚点与 `*` 的键级（2 → ylidene，3 → ylidyne） | `_radical_free_order` → L5 链引擎 |
| `heterane_z` | `int` | 杂原子烃中心元素序数 | `_mononuclear_radical` / `_heterane_parent` → L5 |
| `lambda_n` | `int`/`None` | 单核杂原子烃非标准键数（P-14.1.3） | 同上 → L5 |
| `lambda_atoms` | `dict[int, int]` | 多核均一杂原子链各原子键数 | `_heterane_parent` → L5 |
| `mol` | `Mol` | Mol 引用 | `pack_parent_stem` → 各层 |
| `stem_en` / `stem_zh` | `str` / `str` | 双语词干；阳离子可带 `oxo` 前缀（`oxoazanium`） | `pack_parent_stem` / `_mononuclear_*` → L5 `_ensure_parent_stem` |
| `stem_bare_en` / `stem_bare_zh` | `str` | 不含 `ane`/`烷` 的裸词干（多核均一杂原子链），供链式词干引擎拼词尾或 `a…-triene` | `_heterane_parent` → L5 `_ensure_bridged_stem`/`_ensure_generated_stem` |
| `double_bond`/`double_bonds`、`triple_bond`/`triple_bonds` | `tuple`/`list[tuple]` | 骨架内 C=C / C≡C 端点对 | `_chain_unsat_fields` → L4 |
| `scaffold_id` | `str` | 骨架标识（保留母体 id、`"benzene"`/`"fused"`、生成式 `"hw_mono"`、桥环 `"bridged"`） | `_scaffold_fields` → L4/L5 |
| `scaffold_identity` / `scaffold_match` | `ScaffoldIdentity`/`tuple` | scaffold 对象（`id` 与 `naming_class` 另有 `bridged`/`heterocycle` 取值）；模板原子映射（供 L4 固定编号/算加氢位） | `_scaffold_fields` → L2/L4/L5 |
| `hydro_atoms` | `frozenset[int]` | 加氢位（P-31.2.2）；保留模板取 `hydrogenated_atoms`，生成式 HW 环取 `hantzsch_widman.hydro_atoms` | `_scaffold_fields` → L4 |
| `fused_tree` | `FusedNode` | 稠环拆解树 | `decompose_fused_system` → L4/L5 |
| `bridged_node` / `bridged_nodes` | `BridgedNode` / `tuple[BridgedNode,...]` | L2 下传的并列最优 P-23 桥环候选；L4 裁决后写回选中者 | `_scaffold_fields` → L4 `bridged_numbering` → L5 |
| `spiro_node` / `spiro_nodes`、`fbs_node` / `fbs_nodes` | `SpiroNode` / `FbsNode` 及候选元组 | P-24 螺环 / 组分式螺环候选（`fused_bridged_spiro` 用 `fbs_*` 键） | `_scaffold_fields` → L4 `spiro_numbering`/`fbs_numbering` → L5 |
| `numbering_scaffold` | `dict` | 含 `scaffold_id` 与 `labels`；标签数与 `chain` 不符时不写入 | `_attach_numbering_scaffold`；L4 稠环路径改写为 `"fused"` |
| `owned_atoms` | `frozenset[int]` | 链原子 ∪ 主基团特征原子；含氧酸中心只取后者 | `finalize_parent_ownership` → L3、覆盖台账 |
| `ring_attach_idx` | `int` | 环附着原子（只读键） | 无写入方 → `namer._remap_candidates` |
| `locant_kind` | `str` | 位次省略的名义母体类：`urea`/`thiourea`/`guanidine`/`carbamic_acid` | L5 `_c1_retained` → `_prefix_for`/`_build_prefix` |
| `subs_consumed` | `bool` | 取代基已并入 C1 保留名，前缀不再另加 | L5 `_c1_retained` → `_prefix_for`/`_build_prefix` |
| `stem_generated` | `bool` | `'a'` 前缀生成式词干的标记（自由基取裸词干、位次恒显式） | L5 `_ensure_generated_stem` → `chain_engine` |
| `component_numbering` | `bool` | 稠合组分「仿母体」内部标志（位次由稠合名整体定，不参与收窄；输入侧键） | L5 `fused_namer` → L4 `numbering_engine` |
| `indicated_h_locants` / `indicated_h_forced` | `list[str]` / `bool` | 指示氢位次（P-58.2.1）；`indicated_h_forced = bool(extra) or is_hw_scaffold(scaffold_id)` | L4 `number()` → L5 |
| `ind_h_carbon_ok` | `bool` | 母体是否为可写碳位指示氢的保留杂环名 | L4 `number()` → L5 `join_hydro_prefix` |
| `hydro_prefix` | `tuple[str, str]` | 加氢前缀 `(en, zh)`；位次表达不出时回退 `("", "")` | L4 `number()` → L5 |
| `hydro_fallback` | `bool` | 加氢位被搬进指示氢且前缀已不提它（须写指示氢） | L4 `number()` → L5 |
| `lambda_locants` | `tuple[tuple[str,int,int],...]` | 整体编号下需标注的骨架原子 `((位次, λ键数\|0, δ双键数\|0), …)`（P-14.1.3/P-25.7.2） | L4 `number()` → L5 `_lambda_prefix` |

**PrincipalExpressionFacts**：`group_class`、`multiplicity`、`relation`（`in_skeleton`/`exocyclic`）、`occurrence_ids`、`characteristic_atoms`、`anchor_atoms`（原锚点）、`attachment_atoms`（骨架内附着原子）、`charge_state`。L4 位次读 `attachment_atoms`。

**BridgedNode / BridgeSegment**（`layer2/bridged_system`，P-23）：`atom_ids`、`numbering`（原子 → 位次）、`main_ring`（主环环序）、`ring_segments`（主环两段长）、`main_bridge`、`independent_bridges`/`dependent_bridges`、`ring`/`n_rings`；派生 `secondary_bridges`、`descriptor`（两段 + 主桥 + 各次级桥长度）、`locant_pairs`（次级桥上标位次对）、`scaffold_identity()`（id 与 naming_class 均为 `bridged`）。`BridgeSegment` = `heads`/`atoms`/`numbers`/`locants`。入口 `decompose_bridged_system`、`try_bridged_scaffold`（`_fusion_naming_applies` 为假时接管）。

**SpiroNode / FbsComponent**（`layer2/spiro_system`，P-24）：`SpiroNode` = `scaffold_id`/`atom_ids`/`free_spiro_atoms`/`components`/`descriptor`/`descriptor_superscripts`/`numbering`/`ring`/`n_rings`；`FbsComponent` = `index`/`kind`（`mono_ring`/`fused_ring`/`bridged_ring`）/`atom_ids`/`ring_indices`/`spiro_atoms`/`base_en`/`base_zh`/`bare_en`/`bare_zh`/`sid`/`node`/`numberings`/`order_extra`/`hydro`（`(en, zh)` 加氢前缀）。`SPIRO_SCAFFOLDS = ("mono_spiro", "fused_bridged_spiro")`。

**生成式杂单环**（`layer2/hantzsch_widman`，P-22.2.2）：`HW_ID = "hw_mono"`、`HW_CLASS = "heterocycle"`；`identity` 对未命中保留模板的 3–10 元孤立含杂环发 `ScaffoldIdentity(HW_ID, HW_CLASS, 1, "hetero")`，词干由 `hw_name_from_cycle` 生成（饱和/mancude 两形态）；`hydro_atoms`/`ring_numbering`/`component_key` 供加氢位、环编号与稠合组分编码（`HW_COMPONENT_PREFIX = "hw:"`）复用。

---

## 4. numbered dict（L4 输出 → L5 输入）

产出 `number(parent, substituents)`：`oriented = {**parent, "chain": orient_numbering(...)}`，再由 `locant_calc._pack(oriented, _with_locants(...))` 打包；`number()` 另在 `parent` 上写回指示氢/加氢/λ 字段。

| 键 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `parent` | `dict` | 定向后的 parent（`chain` 为编号顺序），追加 `hydro_fallback`/`lambda_locants`/`indicated_h_locants`/`indicated_h_forced`/`ind_h_carbon_ok`/`hydro_prefix`；L5 可再写 `stem_en`/`stem_zh`/`stem_bare_*`/`stem_generated`/`locant_kind`/`subs_consumed` | `_pack` + `number()` + L5 → L5 |
| `substituents` | `list[dict]` | 参与编号的 subst dict，已注入 `locant`；L5 前缀围栏会就地改写 `en`/`arm_fenced` | `_with_locants` → L5 |
| `fg_locants` | `list[dict]` | `[{kind, locants, omit}]`，只含实际存在位次的 FG | `_fg_locants`（表 `_FG_LOCANTS` 由 `FG_SPECS` 投影）→ L5 |
| `ene_locants` / `yne_locants` | `list[int]` / `None` | 双键 / 三键端点较小侧位次排序列表（跨位次键写作 `x(y)`） | `_unsat_locants` → L5 |
| `omit_ene_locant` / `omit_yne_locant` | `bool` | 烯 / 炔位次省略标志（`locant_calc.omit_unsat`） | 同上 |
| `bridge_self_enclosed` | `bool` | 桥复合前缀名已自含围栏标记位 | L5 `_mononuclear_radical_names` → `namer._chain_meta`、L3 免二次加括号 |

位次来源 `principal_expression_facts.attachment_atoms`（`_locants_for`）；`omit` 由 `locant_calc.omit_fg_locant` 判定，经 `_FG_GROUP`（键同 `FG_SPECS.fg`，值 `alcohol`/`thiol`/`amine`/`ketone`）映射到 principal 类别。`omit_unsat`/`omit_fg_locant`/`suffix_locant_set`/`prefix_locant_set` 均在本模块（`layer4/locant_calc.py`）。

输入侧附加键：`substituents[].fusion`（`bool`，稠合点标记；`orient_numbering`/`_fixed_numbering` 据此让稠合侧落靠前字母，并由 `fused_component_numbering` 生成）；`parent.component_numbering`（稠合组分仿母体内部标志）。对外接口 `numbering_engine.resolve_numbering`（桥环/螺环共用裁决外壳）、`indicated_hydrogen.indicated_hydrogen_atoms`、`indicated_hydrogen.is_hw`。

---

## 5. L3 结构

**ClaimedBlock**：`slot`/`attach_parent`（母体侧附着点）/`root`（取代基侧根原子）/`atoms`。**SideSlot** 四值：`CHAIN_C`、`RING_C`、`AMINE_N`（非芳香非环员电中性 N）、`OTHER`。产出 `iter_claims(mol, owned_atoms)`，按 `(attach_parent, root, slot.value)` 排序；`_try_claim` 丢弃无连接边、或含双键连 owned 非碳重原子的氧的组分（环内 S/P 的 =O 例外，须按 oxo 前缀 claim）。

**CLAIM_KIND** = `{"amine_n": "n_block", "ring_c": "alkyl", "chain_c": "alkyl"}`，兜底 `"side"`（`_claim_kind`）；`NAME_KIND` 登记名称暗含的非烷基 kind（halo 四名）；`sub_from_named` 中附着原子为环员**或**母体链成员且 kind ∈ `N_PREFIX_KINDS` 时改判 `ring_c`（改走位次而非 N- 前缀）。

**SubstituentName**：`claim` / `en` / `zh` / `requires_parentheses`。`SubstituentNamer.name(mol, claim)` 按 `self._backends`（缺省 `[RetainedBackend(), RecursiveBackend(...)]`）取首个命中；`_named` 封装 `(en, zh, paren)`，en 或 zh 为空即失败。

**CoverageLedger**：`owned_atoms` / `named_claims` / `gap` / `overlap`，`complete = not gap and not overlap`；`build_coverage_ledger(mol, *, owned_atoms, names)` 中 `gap = 重原子全集 − owned_atoms`、`overlap` 取归属计数 `> 1` 者。

---

## 6. subst dict（L3 输出 → L4/L5 输入）

产出 `extract_substituents(info, parent, *, cache=None)`；装配函数 `sub_from_named(named, mol, parent=None)`。

| 键 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `kind` | `str` | `NAME_KIND` 优先，否则 `CLAIM_KIND` 映射，兜底 `"side"`；L5 可就地改写为 `"side"`（脲/硫脲/胍/氨基甲酸的 N-取代基改走数字位次通道）或同位素类（`isotope`/`isotope_halo`/`isotope_group`） | `sub_from_named` → `namer._subs_for_numbering`、L4/L5 |
| `en` / `zh` | `str` | 双语取代基名 | 同上 → L4 字母序平局、L5 前缀 |
| `attach_idx` | `int` | 附着于母体的原子索引（= `claim.attach_parent`） | 同上 → L4 位次与编号方向、L5 |
| `atoms` / `n_carbons` | `list[int]` / `int` | 取代基占用的原子索引（升序）/ 取代基碳数 | 同上 → L2 覆盖判定、L5 |
| `paren` | `bool` | 是否需要括号包裹（= `requires_parentheses`） | 同上 → L5 前缀与臂围栏 |
| `o_side` | `bool`（可选） | O 侧烷基臂（parent kind ∈ `ESTER_O_SIDE_KINDS` 或带 `o_idx`，附着原子 O/S） | `_append_named` → `namer._remap_attach`（跳过链重映射）、L5 |
| `locant` | `int` / `str` | 母体骨架位次（保留稠环字母位如 `"4a"`）；不在链上取 0 | L4 `_with_locants` → L5 |
| `arm_fenced` | `bool`（可选） | O-侧臂名已完成围栏的防重入标记 | L5 `_fence_o_side_arms` |
| `n_prime` | `int`（可选） | 酰肼两端 N 的撇号数：`1` = 远端 N（N′，`attach_idx == hydrazide_n_idx`），`0` = 近端 N | L3 `_mark_hydrazide_primes` → L5 `_sub_prime` |
| `fusion` | `bool`（可选） | 稠合点标记（输入侧键，非 L3 产出） | L4 `orient_numbering`/`_fixed_numbering` |

`namer._subs_for_numbering` 只把 `o_side`、`attach_idx ∈ chain`、`kind ∈ N_PREFIX_KINDS` 三者之一的取代基送入 `number()`；其余环外基团留在 numbered `substituents` 中供 L5 拼前缀。

---

## 7. NameResult

`types.NameResult`：`en`、`zh`、`success: bool`、`source = "iupac"`、`time_ms: float`、`meta: dict[str, Any]`。成功由 `assembler.assemble`（8 步：`_ensure_parent_stem` 门控 → `_names_for` 取名 → `join_hydro_prefix` → `join_ring_cation_suffix` → `join_isotope_descriptor` → `_ring_oxide_suffix` → `join_kind_name(kind, _prefix_for(...), names, numbered)` → `_chalcogen_nitrile_alias` → `join_anion_names` → `join_ez_prefix` → `join_rs_prefix`）实例化、`namer._ok_result` 补链元数据；失败由 `namer._fail`/`assembler._fail`。

| `meta` 键 | 类型 | 含义 | 写入方 |
|---|---|---|---|
| `reason` / `n_carbons` / `kind` | `str` / `int` / `str` | 失败原因（parse / no_assemblable_candidate / unsupported / simple_species）；`unsupported` 另带 `n_carbons`/`kind` | `namer._fail`、`assembler._unsupported`、`namer._name_mol` |
| `parent_chain` | `list[int]` | 母体链原子索引；写缓存前经 `_canonical_result` 规范排序 | `namer._chain_meta` |
| `parent_kind` | `str` | 母体 `kind`；`∈ OXO_CENTER_KINDS` 时跳过盐后缀 | `namer._chain_meta` → `_apply_salt_suffix` |
| `parent_labels` / `bridge_self_enclosed` | `list` / `bool` | 整体编号标签（稠环桥头 `3a`/`6a`，长度与 `chain` 不符时为空表）；桥前缀已自含围栏标记 | `namer._chain_meta` |
| `parent_substituent_count` | `int` | `len(numbered["substituents"])`，供 L3 判词干是否复合 | `namer._ok_result` |
| `p44_1_1_key` / `p45_2_2_key` | `tuple` | P-44.1.1 后缀 / P-45.2.2 前缀位次集合（`locant_calc.suffix_locant_set`/`prefix_locant_set`） | `namer._assemble_candidate`；`_best_hit` 取最小者 |
| `fallback` | `str` | 恒为 `"no_coverage_gate"` | `namer._try_phase` |
| `salt` | `dict` | L0 盐元数据（非空且成功时写入） | `namer._name_mol` |

顶层入口 `SMILESNNamer.name(smiles)`：查 `CommonNameCache`，未命中则 `memo.begin_run()` → `_pipeline` → `_name_mol`，成功结果经 `_canonical_result` 写回缓存。`_name_mol` 顺序为：`_simple_species` 查表短路 → `dissociate_salt` → `analyze` → 注入 `root_ctx`/`salt` → `_run_candidates` → `_apply_salt_suffix` → 多片段/失败回落 `_name_components`（`pair_unique_ions` → `_merge_repeats`）。

---

相关页面：[[architecture/layer0-preprocessor]]、[[architecture/layer1-analyzer]]、[[architecture/layer2-parent-selector]]、[[architecture/layer3-substituents]]、[[architecture/layer4-numbering]]、[[architecture/layer5-name-assembly]]、[[concepts/atom-ownership]]、[[concepts/bilingual-naming]]、[[concepts/functional-group-priority]]。
