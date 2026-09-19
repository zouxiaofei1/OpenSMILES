# 核心数据合约 (Core Data Contracts)

> NamePredict 六层命名流水线跨层数据结构的字段规格；字段名与语义以 `src/namepredict` 源码为唯一依据。

## 端到端数据流图

```text
SMILES ─preprocess─► Mol ─dissociate_salt─► (organic Mol, salt)
  │ analyze(organic)
  ▼
info dict ─namer._name_mol 注入 root_ctx/salt─►
  │ select_parent(info)
  ▼
list[parent dict] ─finalize_parent_ownership─► owned_atoms
  │ extract_substituents(info, parent)
  │   ├ iter_claims → ClaimedBlock       (layer3/claimable_block)
  │   └ SubstituentNamer.name → SubstituentName
  ▼
list[subst dict] ─number(parent, subs)─► numbered dict
  │   └ 桥环：bridged_nodes(L2 候选) → bridged_node(L4 写回) → L5 描述符/上标
  │ assemble(numbered)
  ▼
NameResult ─namer._ok_result─► name(smiles)
```

覆盖门控 `build_coverage_ledger(..., names=[])` 只判 gap（`namer._prepare_candidate`）。

---

## 1. info dict（L1 输出）

产出 `analyzer.analyze(mol)`（基础 + `_collect_fgs` + `_ring_meta`），`namer._name_mol` 注入 `root_ctx`、`salt`；清单取用入口 `inventory_from_info`（缺失抛 `KeyError`）。

| 字段 | 类型 | 含义 | 生产 → 消费 |
|---|---|---|---|
| `mol` | `Mol` | 去盐后的有机部分 | `_info` → 各层 |
| `carbon_ids` / `n_carbons` | `list[int]` / `int` | 碳原子索引及数量 | `_info` → L1/L2 骨架枚举 |
| `double_bonds` / `triple_bonds` | `list[dict]` | `{"c1","c2"}`，C=C / C≡C（非芳香） | `_filter_bond_entries` → L2 `_chain_polys` |
| `fg_inventory` | `FunctionalGroupInventory` | 官能团事实的唯一出口 | `build_inventory` → L2/L3 |
| `rings` / `n_rings` / `has_ring` | `list[dict]` / `int` / `bool` | SSSR 环及数量、是否含环 | `_ring_meta` → `namer._prepare_candidate` |
| `ring_systems` / `n_ring_systems` | `list[dict]` / `int` | 稠合 + 螺环合并的环系（`atom_ids`/`sssr_indices`/`fusion_edges` 等） | `build_ring_systems` → L2 `_scaffold_fields`、L4 |
| `root_ctx` | `tuple[Mol, list[int]]` | `(根分子, 本分子原子→根索引)`；顶层为 `(organic, range(n))` | `namer._name_mol` → L3 在根分子上重算 R/S |
| `salt` | `dict` | L0 盐元数据（`metal`/`metal_zh`/`n_metal`/…） | `namer._name_mol` → L2 盐门控（写 `salt_meta`）、盐后缀 |

---

## 2. FG occurrence / inventory

`FunctionalGroupOccurrence`（`layer1/functional_group_inventory.py`）：

| 字段 | 类型 | 含义 |
|---|---|---|
| `id` | `str` | `"{fg}:{i}"`，如 `"acid:0"` |
| `group_class` | `FunctionalGroupClass` | 14 个注册 FG 类别 + `NONE = "alkane"` |
| `characteristic_atoms` | `frozenset[int]` | 默认 `center_idx` ∪ `surr_idx` |
| `parent_anchors` | `frozenset[int]` | 由 `FgSpec.anchors` 声明的 key 从 payload 取出 |
| `payload` | `dict` | L1 原始条目 |
| `demoted` | `bool` | P-41 降级前缀叶（carboxy/cyano），碳排除出主链 |

`FunctionalGroupInventory`：字段 `entries`；方法 `occurrences(group_class)`（该类未降级条目）、`demoted_entries()`；构建入口 `build_inventory(lists, mol, demoted)`。

**P-41 仲裁**（`_arbitrate_parts`，先于清单构建）：`_SUPPRESSIBLE`（acid/ester/acyl_halide/amide/nitrile/aldehyde）遇 `p41` 更小者整组退出；`_LEAF_DEMOTED`（acid/nitrile）标 `demoted`、碳排除出主链，其余组合 FG 整组清空；`_PRESENCE_SKIP`（oxoacid/sulfonamide）不判存在性。返回 `(parts, demoted_ids)`。

**payload 键**：

| 类别 | payload 键 |
|---|---|
| `radical` | `center_idx`、`surr_idx = []` |
| `acyl`、`acid`、`ester`、`acyl_halide`、`amide`、`nitrile`、`aldehyde`、`ketone`、`alcohol`、`thiol`、`amine` | `center_idx`、`surr_idx` |
| `oxoacid`、`sulfonamide` | `oxo_z`、`center_idx`、`oxo_kind`、`n_oh`、`n_om`；**仅碳锚定条目**另有 `surr_idx`，只收阴离子氧（中心自任母体时无该键） |

- `n_oh` = 中性含 H 氧数；`n_om` = 阴离子氧数；`n_arms` 不在 payload，L2 取 `payload.get("n_arms", 0)`。
- `oxo_kind`：P — `phosphate`/`phosphonate`；S — `sulfate`/`sulfonic`/`sulfonate`/`sulfonamide`/`sulfonyl_chloride`；`_OXO_CLASS_BY_KIND` 只分走 `sulfonamide`。
- 锚点 key：`center_idx` 用于除 `alcohol`/`thiol`/`amine`（取 `surr_idx`）外各类；特征原子例外 `FG_ATOM_FNS`（oxoacid/sulfonamide）取中心 P/S + 非碳邻居 + 锚点碳。

---

## 3. parent dict（L2 输出）

产出 `select_parent(info)` → `list[dict]`：P-44 评分经 P-45.2 / P-45.2.1 重排的并列最优组，整组逐个尝试 L4–L5。构建链 `rule_driven_parent_candidates` → `express_*_principal` → `_parent_dict` → `pack_parent_stem` → `finalize_parent_ownership`。

| 字段 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `kind` | `str` | FG 类别名 / `oxo_kind` / `"radical"` / `"acyl"` / `"alkane"` / `"bridged"` | `_parent_dict` → L4/L5 |
| `chain` / `n_carbons` | `list[int]` / `int` | 骨架原子顺序及长度 | `_parent_dict`；L4 `number()` 覆盖为编号顺序 |
| `covered_principal_ids` | `tuple[str,...]` | 本骨架表达的主基团 id | `_parent_dict` → L2 |
| `principal_occurrences` | `tuple[FunctionalGroupOccurrence,...]` | 主基团 occurrence 全集 | `_parent_dict` → L2 |
| `principal_group_count` | `int` | 本骨架覆盖数 | `_parent_dict` → L2 |
| `principal_expression_facts` | `PrincipalExpressionFacts` | 主基团表达事实（见下） | `_facts` → L3/L4 |
| `radical_c_idx` / `acyl_c_idx` | `int` | 固定 locant 1 单锚点字段（仅 RADICAL/ACYL） | `_semantic_anchor_fields` → L4 |
| `anion` | `bool` | 酸全阴离子（仅 ACID/OXOACID） | `_expression_flags` → L5 |
| `oxo_kind` | `str` | 含氧酸 kind | `_chain_oxoacid_fields` → L5 |
| `n_oh` / `n_om` / `n_arms` | `int` | 酸式 H 氧 / 阴离子氧 / O–R 臂数 | 同上 → L5 |
| `salt_meta` | `dict`/`None` | 盐门控后的 L0 盐元数据；`n_om > 0` 须 `n_metal == n_om` | 同上 → L5 |
| `o_idx` / `alkoxy_n` | `int` | 酯侧杂原子索引（硫酯为 S）/ 臂计数 | `_chain_ester_fields` → L3/L5 |
| `thio_side` / `hal_idx` / `hal_z` | `bool`/`int`/`int` | 硫代酯标记；酰卤卤原子索引与原子序 | `_thio_fields` / `_chain_acyl_halide_fields` → L5 |
| `radical_anchor_element` / `radical_ylidene` | `str`/`bool` | 杂原子锚点元素 / 碳锚点双键自由价 | `_mononuclear_radical` / `_radical_ylidene` → L5 |
| `mol` / `stem_en` / `stem_zh` | `Mol`/`str`/`str` | Mol 引用；双语词干 | `pack_parent_stem`、`_mononuclear_radical`、L5 `_ensure_parent_stem` |
| `stem_bare_en` / `stem_bare_zh` | `str` | 不含 `ane`/`烷` 的裸词干，供链式词干引擎拼词尾或 `a…-triene` | L5 `_ensure_bridged_stem`/`_ensure_generated_stem` → `_names_for` |
| `bridged_nodes` / `bridged_node` | `tuple[BridgedNode,...]` / `BridgedNode` | L2 下传的并列最优 P-23 桥环候选；L4 裁决后写回选中者 | `_scaffold_fields` → L4 `bridged_numbering` → L5 |
| `subs_consumed` / `locant_kind` | `bool` / `str` | 取代基已并入 C1 保留名，前缀不再另加 / 位次省略的名义母体类：`urea`/`thiourea`/`guanidine`/`carbamic_acid` | L5 `_c1_retained`/`_names_for` → `_prefix_for`/`_build_prefix` |
| `double_bond`/`double_bonds`、`triple_bond`/`triple_bonds` | `tuple`/`list[tuple]` | 骨架内 C=C / C≡C 端点对 | `_chain_unsat_fields` → L4 |
| `scaffold_id` | `str` | 骨架标识（保留母体 id 或 `"benzene"`/`"fused"` 等） | `_scaffold_fields` → L4/L5 |
| `scaffold_identity` / `scaffold_match` / `typed_ring_expression_supported` | `ScaffoldIdentity`/`tuple`/`bool` | scaffold 对象（`id` 与 `naming_class` 另有 `bridged` 取值）；模板原子映射；能否 typed 表达 | `_scaffold_fields` → L2/L4/L5 |
| `hydro_atoms` | `frozenset[int]` | 加氢位（P-31.2.2） | `_scaffold_fields` → L4 |
| `fused_tree` | `FusedNode` | 稠环拆解树 | `decompose_fused_system` → L4/L5 |
| `numbering_scaffold` | `dict` | 含 `scaffold_id` 与 `labels`；标签数与 `chain` 不符时不写入 | `_attach_numbering_scaffold`；L4 稠环路径改写为 `"fused"` |
| `owned_atoms` | `frozenset[int]` | 链原子 ∪ 主基团特征原子；含氧酸中心只取后者 | `finalize_parent_ownership` → L3、覆盖台账 |
| `ring_attach_idx` | `int` | 环附着原子（只读键） | 无写入方 → `namer._remap_candidates` |
| `indicated_h_locants` / `indicated_h_forced` | `list[str]` / `bool` | 指示氢位次（P-58.2.1）；来自保留名未隐含位 | L4 `number()` → L5 |
| `hydro_prefix` | `tuple[str, str]` | 加氢前缀 `(en, zh)`；位次表达不出时回退 `("", "")` | L4 `number()` → L5 |

**PrincipalExpressionFacts**：`group_class`、`multiplicity`、`relation`（`in_skeleton`/`exocyclic`）、`occurrence_ids`、`characteristic_atoms`、`anchor_atoms`（原锚点）、`attachment_atoms`（骨架内附着原子）、`charge_state`。L4 位次读 `attachment_atoms`。

**BridgedNode / BridgeSegment**（`layer2/bridged_system`，P-23）：`atom_ids`、`numbering`（原子 → 位次）、`main_ring`（主环环序）、`ring_segments`（主环两段长）、`main_bridge`、`independent_bridges`/`dependent_bridges`、`ring`/`n_rings`；派生 `secondary_bridges`、`descriptor`（两段 + 主桥 + 各次级桥长度）、`locant_pairs`（次级桥上标位次对）、`scaffold_identity()`（id 与 naming_class 均为 `bridged`）。`BridgeSegment` = `heads`/`atoms`/`numbers`/`locants`。入口 `decompose_bridged_system`、`try_bridged_scaffold`（稠合不适用、丢环或 P25 产名失败时接管）。

---

## 4. numbered dict（L4 输出 → L5 输入）

产出 `number(parent, substituents)`：`oriented = {**parent, "chain": orient_numbering(...)}`，再由 `locant_calc._pack(oriented, _with_locants(...))` 打包。

| 键 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `parent` | `dict` | 定向后的 parent（`chain` 为编号顺序），追加 `indicated_h_locants`/`indicated_h_forced`/`hydro_prefix`；L5 可再写 `stem_en`/`stem_zh`/`stem_bare_*` | `_pack` + `number()` → L5 |
| `substituents` | `list[dict]` | 参与编号的 subst dict，已注入 `locant`；L5 前缀围栏会就地改写 `en`/`arm_fenced` | `_with_locants` → L5 |
| `fg_locants` | `list[dict]` | `[{kind, locants, omit}]`，只含实际存在位次的 FG | `_fg_locants`（表 `_FG_LOCANTS` 由 `FG_SPECS` 投影）→ L5 |
| `ene_locants` / `yne_locants` | `list[int]` / `None` | 双键 / 三键端点较小侧位次排序列表 | `_unsat_locants` → L5 |
| `omit_ene_locant` / `omit_yne_locant` | `bool` | 烯 / 炔位次省略标志（`omit_locants.omit_unsat`） | 同上 |
| `bridge_self_enclosed` | `bool` | 桥复合前缀名已自含围栏标记位 | L5 `_mononuclear_radical_names` → `namer._chain_meta`、L3 免二次加括号 |

位次来源 `principal_expression_facts.attachment_atoms`（`_locants_for`）；`omit` 由 `omit_locants.omit_fg_locant` 判定，经 `_FG_GROUP`（键同 `FG_SPECS.fg`，值 `alcohol`/`thiol`/`amine`/`ketone`）映射到 principal 类别。

---

## 5. L3 结构

**ClaimedBlock**：`slot`/`attach_parent`（母体侧附着点）/`root`（取代基侧根原子）/`atoms`。**SideSlot** 四值：`CHAIN_C`、`RING_C`、`AMINE_N`（非芳香非环员 N）、`OTHER`。产出 `iter_claims(mol, owned_atoms)`，按 `(attach_parent, root, slot.value)` 排序；`_try_claim` 丢弃无连接边、或含双键连 owned 非碳重原子的氧的组分。

**CLAIM_KIND** = `{"amine_n": "n_block", "ring_c": "alkyl", "chain_c": "alkyl"}`，兜底 `"side"`（`_claim_kind`）；`NAME_KIND` 登记名称暗含的非烷基 kind（halo 四名）；附着原子为环员且 kind ∈ `N_PREFIX_KINDS` 时改走环上位次。

**SubstituentName**：`claim` / `en` / `zh` / `requires_parentheses`。`SubstituentNamer.name(mol, claim)` 按 `self._backends`（缺省 `[RetainedBackend(), RecursiveBackend(...)]`）取首个命中；`_named` 封装 `(en, zh, paren)`，en 或 zh 为空即失败。

**CoverageLedger**：`owned_atoms` / `named_claims` / `gap` / `overlap`，`complete = not gap and not overlap`；`build_coverage_ledger(mol, *, owned_atoms, names)` 中 `gap = 重原子全集 − owned_atoms`、`overlap` 取归属计数 `> 1` 者。

---

## 6. subst dict（L3 输出 → L4/L5 输入）

产出 `extract_substituents(info, parent, *, cache=None)`；装配函数 `sub_from_named(named, mol)`。

| 键 | 类型 | 含义 | 写入方 → 读取方 |
|---|---|---|---|
| `kind` | `str` | `NAME_KIND` 优先，否则 `CLAIM_KIND` 映射，兜底 `"side"`；L5 可就地改写为 `"side"`（脲/硫脲/胍/氨基甲酸的 N-取代基改走数字位次通道） | `sub_from_named` → `namer._subs_for_numbering`、L4/L5 |
| `en` / `zh` | `str` | 双语取代基名 | 同上 → L4 字母序平局、L5 前缀 |
| `attach_idx` | `int` | 附着于母体的原子索引（= `claim.attach_parent`） | 同上 → L4 位次与编号方向、L5 |
| `atoms` / `n_carbons` | `list[int]` / `int` | 取代基占用的原子索引（升序）/ 取代基碳数 | 同上 → L2 覆盖判定、L5 |
| `paren` | `bool` | 是否需要括号包裹（= `requires_parentheses`） | 同上 → L5 前缀与臂围栏 |
| `o_side` | `bool`（可选） | O 侧烷基臂（parent kind ∈ `ESTER_O_SIDE_KINDS` 或带 `o_idx`，附着原子 O/S） | `_append_named` → `namer._remap_attach`（跳过链重映射）、L5 |
| `locant` | `int` / `str` | 母体骨架位次（保留稠环字母位如 `"4a"`）；不在链上取 0 | L4 `_with_locants` → L5 |
| `arm_fenced` | `bool`（可选） | O-侧臂名已完成围栏的防重入标记 | L5 `_fence_o_side_arms` |

`namer._subs_for_numbering` 只把 `o_side`、`attach_idx ∈ chain`、`kind ∈ N_PREFIX_KINDS` 三者之一的取代基送入 `number()`；其余环外基团留在 numbered `substituents` 中供 L5 拼前缀。

---

## 7. NameResult

`types.NameResult`：`en`、`zh`、`success: bool`、`source = "iupac"`、`time_ms: float`、`meta: dict[str, Any]`。成功由 `assembler._ok`（`assemble` 入口）实例化、`namer._ok_result` 补链元数据；失败由 `namer._fail`/`assembler._fail`。

| `meta` 键 | 类型 | 含义 | 写入方 |
|---|---|---|---|
| `reason` / `n_carbons` / `kind` | `str` / `int` / `str` | 失败原因（parse / no_assemblable_candidate / unsupported）；`unsupported` 另带 `n_carbons`/`kind` | `namer._fail`、`assembler._unsupported` |
| `parent_chain` | `list[int]` | 母体链原子索引；写缓存前经 `_canonical_result` 规范排序 | `namer._chain_meta` |
| `parent_kind` | `str` | 母体 `kind`；`∈ OXO_CENTER_KINDS` 时跳过盐后缀 | `namer._chain_meta` → `_apply_salt_suffix` |
| `parent_labels` / `bridge_self_enclosed` | `list` / `bool` | 整体编号标签（稠环桥头 `3a`/`6a`，长度与 `chain` 不符时为空表）；桥前缀已自含围栏标记 | `namer._chain_meta` |
| `parent_substituent_count` | `int` | `len(numbered["substituents"])`，供 L3 判词干是否复合 | `namer._ok_result` |
| `p44_1_1_key` / `p45_2_2_key` | `tuple` | P-44.1.1 后缀 / P-45.2.2 前缀位次集合（`candidate_keys`） | `namer._assemble_candidate`；`_best_hit` 取最小者 |
| `fallback` | `str` | 恒为 `"no_coverage_gate"` | `namer._try_phase` |
| `salt` | `dict` | L0 盐元数据（非空且成功时写入） | `namer._name_mol` |

顶层入口 `SMILESNNamer.name(smiles)`：查 `CommonNameCache`，未命中则 `memo.begin_run()` → `_pipeline` → `_name_mol`，成功结果经 `_canonical_result` 写回缓存。
