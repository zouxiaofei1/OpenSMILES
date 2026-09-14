# layer0-4 代码分析与精简报告

范围：`src/namepredict/layer0` ~ `layer4`（约 6200 行）+ `constants.py` / `types.py` / `namer.py`。
结论：代码行数是负资产，本报告的删减项合计约 **700 行**（占比约 12%），另有 3 处真实缺陷与 10 处高收益热点。

> **已执行**：§1.1 酸酐路径与 §2 首行桥环路径已删除（`ring_systems.py` 301 → 209 行，全仓 anhydride / bridge 标识符清零）。
> benchmark 实测 EN 58.5% / ZH 76.3% / dual 58.2%，`IMPROVE=1 REGRESS=0`。
> 原先崩溃的十二酸酐 `CCCCCCCCCCCC(=O)OC(=O)CCCCCCCCCCC` 由 fail → ok，预测 `dodecanoyl dodecanoate` 与 gold 一致；
> 其余酸酐分子现按酯命名（`acetyl acetate` / `propanoyl propanoate`），不再抛 `NameError`。

---

## 一、真实缺陷

### 1.1 酸酐路径是活的崩溃点（已实测复现）

`CC(=O)OC(C)=O`、`CCC(=O)OC(=O)CC` 走 `analyze()` 直接抛 `NameError: name 'other' is not defined`。

- `analyzer.py:301` — `c1, c2 = sorted((atom.GetIdx(), True))` 之后引用了不存在的 `other`
- `analyzer.py:313-321` — `_anhydride_entries` 构造 `e` 后从不 `append`，即使修好 NameError 也永远检不出酸酐
- `analyzer.py:147-152` `_anhydride_o_of`、`:179-188` `_is_anhydride_carbon`、`:305-311` `_ald_blocked` 内分支 —— 全部只服务这条已崩的路径

### 1.2 `SideSlot.AMIDE_N` 全链路不可达

- `claimable_block.py:29-33` — `_is_amide_n` 在元素检查后无条件 `return False`，实现被挖空，docstring 仍描述原规则
- 连带 `constants.py:154` `CLAIM_KIND["amide_n"]` 成为死登记，`derive_slot` 首个分支恒不触发

### 1.3 `RING_TEMPLATES` 越界

- `ring_geometry.py:19-21` — 只建 `range(3, 20)`
- `fused_orientation.py:148,225` — 裸下标 `RING_TEMPLATES[n]`，≥20 元大环直接 KeyError

### 1.4 多片段分子漏做质子重定位

- `charge.py:93-98` — donor 只取排序后的 `donors[0]`，其所在片段无弱受体即 `break`，不再尝试同片段有受体的次优先 donor。多片段分子漏重定位。

### 1.5 P-44 候选级评分实际不存在

- `parent_selector.py:21` `principal_key` 无任何调用者
- `parent_selector.py:25-27` `_rank_candidates` 是 `key=lambda c: []` 的恒等排序（no-op），docstring 却称"按评分降序"
- `P44Facts` / `ParentCandidate` 在 src 零消费；`P44Facts.principal_group_class: int` 实际传 `None`

候选顺序目前完全由收集顺序决定。

---

## 二、可删除（零消费者，纯减行）

| 位置 | 内容 | 量 |
|---|---|---|
| `ring_systems.py:83-193` | 整套桥环机制 10 个函数（`_ring_adjacent` / `_non_adjacent_pairs` / `_bridgeheads` / `_walk_path` / `_paths_between` / `_dedup_paths` / `_bridge_paths` / `_bridge_info` / `_try_bridge_component` / `_compute_bridged_info`），产出的 `bridgeheads` / `bridge_lengths` / `bridge_paths` 全仓零消费者 | ~110 行 |
| `ring_systems.py:74-81,213-243` | `_hetero_atoms` 及合并逻辑；`spiro_in` / `is_bridged` / `mol` 三个死参数；恒 `None` 的 `topology` / `is_aromatic_mancude`（`server/backend/locants_svg.py:36` 读 `is_aromatic_mancude` 因此恒假） | ~30 行 |
| layer2 死字段 | `ScaffoldSpec.sub_rules` / `principal_slots`；`NumberingPolicy.anchors` / `substitutable` + 恒不可达分支；`next_rule` / `unsupported_ids`；`_candidate_key`；`relative_stereo`；`RetainedEntry`；`ScaffoldId` | ~70 行 |
| layer3 转发壳 | `_yl_from_sub`、`_ordered`、`_retained_name`、`_carbon_slot`、`namer._claim_from_sub`、`claim_extract` 的 `existing=` / `_covered_atoms` / `covered` 整条路径 | ~45 行 |
| layer4 死参数/分支 | `FIXED_START_KEYS` 未用导入；`anchor_as_principal` 死分支；`has_ene` / `has_yne` 恒 None 参数；`preferred_orientations` 未用的 `mol` 形参；`ring_geometry` 的 `right_edge_vertical=False` 分支 | ~30 行 |
| layer0 / constants | `OH_KINDS` / `AMINE_KINDS` 零引用；`salt.py:7` 未用 `O` 导入；`charge.py:75` 注释残留；`charge.py` 的 `min_oxo` 恒 1 形参；`ACID_KIND_PRIO`（可被 `DONOR_KIND.index` 替代） | ~20 行 |
| layer1 | `fg_registry` 的 `oh_parent` / `nh2_parent` / `oxo_parent` 零读取；`_heavy(mol, a)` 未用形参；`_is_acyl_head` 恒真条件 | ~20 行 |

---

## 三、可合并：同一件事的多份实现

### 位次映射共 10 份（最值得统一）

`locant_calc.py:24`、`locant_calc.py:59`、`candidate_keys.py:18`、`numbering.py:34/59/75`、`indicated_hydrogen.py:45`、`numbering_engine.py:239`、`fused_numbering.py:141`、`layer5/stereo.py:111`、`layer5/assembler.py:508`

全部是 `labels[chain.index(a)] if len(labels) == len(chain) else index + 1`。L5 的 docstring 自认"与 L4 `_atom_locant` 同约定"，说明口径已漂移。收敛到 `locant_calc._atom_locant`，不新增依赖方向。

### 候选收窄骨架 11 份

`numbering_engine._narrow:71`、`fused_numbering._keep:161`、`fused_system._keep_best:79` / `_gj:127`、`parent_skeleton` 七个 `keep_*` —— 同一"取 max/min key 保并列"，三种 None 处理各写一遍。

### 其余高频重复

- **fusion_edges 子环集重映射** — `fused_system.py:159-161`、`numbering_engine.py:278-280`、`fused_namer.py:29-31` 逐字重复 → 抽 `sub_system_edges(rset, fusion_edges)`
- **环覆盖计数** — `parent_skeleton.py:179`、`principal_expression.py:157`、`ring_scaffold.py:470`、`numbering_engine.py:263` 四处一行不差 → `ring_count_within(mol, atom_ids)` 放进 `ring_systems`
- **sanitize + AssignStereochemistry 尾巴** — `charge.py:61-65`、`tautomer.py:65-69`、`preprocessor.py:19-20`；`numbering_engine.py:100` 的 flags 已漂移成 `cleanIt=True` → 抽 layer0 私有 helper
- **酯烷氧基谓词 4 个 → 1 个** — `analyzer.py:155-165` 与 `acyl_halide.py:15-25` 各写一份同名包装，归一后行为等价
- **`_sssr` / `sssr_rings`** — `ring_systems.py:10-16` 同一实现两个名字
- **`_components` / `_group_by_root`** — `ring_systems.py:64-72` vs `:265-270` 同一算法两份
- **`_unsat_bonds` 同名两义** — `numbering_engine.py:145`（parent dict → 键对）vs `locant_calc.py:70`（oriented dict + "ene"/"yne" → 列表）
- **位次元组三种形态** — `numbering_engine._locant_set`（tuple | None）、`fused_numbering._locant_tuples`（tuple）、`fused_system._locant_tup`
- **标注与实现不符** — `locant_calc.py:34` 标注 `-> list[int]`，稠环实际返回 `"4a"` 字符串，已扩散到全部消费点

---

## 四、可优化：性能热点

按收益排序，均不改返回值。

1. **`preferred_orientations` 无记忆化** — `numbering_engine.py:281,388`。L2 `fused_system._numbered_locants:153` 对每个母体候选各调一次（几何摆放 + 外周编号），是整条管线最大热点。用现成 `memo.by_key`。
2. **`build_ring_systems` 无记忆化** — `ring_systems.py:289`。一次命名被完整重算 2 次（`analyzer._ring_meta` + L2/L4/L5）。加 `memo.by_mol("ring_systems", ...)`。
3. **`_elem_sig` 无记忆** — `ring_scaffold.py:332-337`。单分子实测被调 **108 次**（模板侧有 `_TEMPLATE_ELEM`，分子侧没有）。
4. **`_match_with_map` 全表扫 83 模板 + 82 氢化模板** — `ring_scaffold.py:410-440`。实测四环分子调用 17 次，fused 拆解占 L2 耗时 60%+。
5. **同一骨架算两遍 match** — `ring_scaffold.py:342` 与 `principal_expression.py:199-210`；`numbering_engine.py:190-201` 又重跑 `GetSubstructMatches`，而 L2 已存入 `parent["scaffold_match"]`。
6. **L2 为每个候选跑完整 L4 编号两次** — `fused_system.py:153-171` 评分时算一遍，`numbering_engine.py:251-311` 胜出后又从零算一遍。
7. **`fused_numbering._keep` 每候选算 `_locant_tuples` 3 次** — `:161-164`；`_alpha_key` / `_rs_locant_key` 各 2 次；`_locant_tuples` 内还是 `chain.index` 的 O(n)/原子。
8. **`indicated_hydrogen.py:14-18` 每次调用 `Chem.Mol(mol)` 全分子复制 + Kekulize** — 一次命名内 `saturated_ring_atoms` 被调 3-4 次。
9. **`claimable_block` 同一 root 的 `cut_block` BFS 跑两遍** — `:81` 与 `:132-133`；且 `_attach_parents_of` / `_canonical_edge` / `_has_dbl_o_edge` 把同一份原子集扫 3 遍。
10. **`analyzer` 羰基谓词重复计算** — 实测 23 原子分子一次 `analyze`：`_is_single_c_oh` 65 次、`_is_carboxylate_o` 57 次、`_has_double_bonded_o` 26 次。5 个谓词对同一原子重算同一事实，可在 `_detect_parts` 内建 per-atom 事实 dict。

---

## 五、可改数据驱动

- **`constants.py:41-42`** — `P25_SENIOR` 与 `P145_SENIOR` 两串 21 元素手工维护、实际只差 `N` 的位置 → `P25_SENIOR = (N,) + tuple(z for z in P145_SENIOR if z != N)`
- **`FG_PARTS_KEY`** — `constants.py:114` 与 `analyzer.py:439` 的派生式实测完全相同只差 `phosphate`。改用派生式，同时消掉"重复登记"与"磷酸不参与 P-41 仲裁"两个问题
- **`locant_calc.py:105` `_FG_GROUP`** — 与同文件 `:140` `_FG_LOCANTS`（已由 `FG_SPECS.locant_kind` 派生）是同一映射两份 → 删 `_FG_GROUP`
- **`kind_registry.py` 的 `_REG`** — 是 `ScaffoldSpec` 的全量复投影（83 条完全一致）→ 删注册表，`pack_parent_stem` 直接 `get_spec(kind)`
- **`layer5/chain_engine._KIND_TABLE`** — 其 `fg` / `mult_ok` 与 `FG_SPECS` 的 `locant_kind` / `multi` 同义 → 由 `FG_SPECS` 派生
- **`ring_geometry.py:25-54`** — 变形环模板硬编码仅 `n∈{5,7}` 且 `exit_idx∈{2,3}`，n=9/11 返回 None → 改 `(n, exit_idx) → 坐标` 数据表
- **`_chain_orient.py`（32 行三层包装）、`candidate_keys.py`（25 行两函数）** — 溶解，同时消掉一处 lazy import 规避的循环依赖
- **`anchored_table.py:29 与 :31`** — 重复登记 `"nitro"` 键（后者静默覆盖）
- **元素字面量 `7/6/8/1` 散落 35+ 处** — `claimable_block`、`analyzer`、`parent_skeleton`、`numbering_engine` 等。`constants.py` 已有 `N/C/O/S` 符号，换符号后「元素支持集」可 grep

---

## 六、文件级合并

13 → 5：

| 合并后 | 组成 | 理由 |
|---|---|---|
| `parent_select.py` | parent_selector + candidates + principal_parent + parent_ownership | 一条无分支单向转发链，没有任何第二种调用路径 |
| `principal_expression.py` | principal + principal_expression + ring_expression_policy | principal 是已被 import 的注册表，policy 是 17 行数据表且只有一个消费者 |
| `ring_scaffold.py` | ring_scaffold + kind_registry | kind_registry 是 ScaffoldSpec 的全量复投影，典型"只转发的注册表" |
| `parent_skeleton.py` | parent_skeleton + chain_walk | chain_walk 的 3 个被引函数只服务 P-44 骨架枚举 |
| `fused_system.py` | 保持独立 | 稠环拆解算法自成一体 |

layer3：`substituent_extractor.py`（整文件是转发壳）并入 `claim_extract.py`；`SubstituentNamer` 的 Protocol + 2 个 Backend + 工厂（约 60 行，`name` 属性生产代码从不读）塌缩为单个函数（约 20 行）。

---

## 七、import 拓扑问题

- **layer2 ↔ layer4 双向环** — `fused_system.py:119,120,162,163` → layer4；`numbering_engine.py` 六处、`numbering.py:117` → layer2.ring_scaffold
- **低层 import 高层** — `layer2/parent_selector.py:33` → `layer3.claimable_block.iter_claims`
- **tools 依赖流水线层** — `tools/anchored_table.py:10` → `layer3.submol_build`，与 `tools/__init__.py` 自称"不依赖任何流水线层"矛盾
- **namer ↔ layer3 环** — `layer3/as_substituent.py:17` import `layer5.stereo` 私有函数，L3 依赖 L5 才有的 R/S 校正
- **模块内惰性 import 20+ 处** — `principal_expression.py` 5 处、`analyzer.py` 3 处、`claim_extract.py:67` 等，多数不存在循环依赖，可直接上提

---

## 八、建议执行顺序

1. **修缺陷** — 酸酐路径（删或修）、`_is_amide_n` 挖空、`RING_TEMPLATES` 越界、`charge.py` 多片段重定位
2. **删死代码** — 约 350 行零消费者代码（桥环 110 行是最大单块）
3. **统一位次映射与环覆盖计数** — 约 -80 行，消除口径漂移
4. **加 3 处 memo** — `preferred_orientations` / `build_ring_systems` / `_elem_sig`，各 3 行，收益最大
5. **layer2 文件合并** — 约 -300 行
