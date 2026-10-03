# Layer4: 编号与位次 (Numbering & Locants)

> 源文件 10 个 `.py`（2,072 行）｜对外接口 `numbering.number`

## 概述

Layer4 承接 Layer2 选定的母体 `parent` 与 Layer3 提取的取代基 `substituents`，按 P-14.4 / P-23.3.2 / P-24 / P-25.3.3 给母体骨架定向编号，并把编号结果物化为 Layer5 可直接消费的位次数据。

职责边界：

- 只决定骨架的「原子顺序 / 位次」，不生成名称文本、不选母体、不改动取代基集合。
- 编号不适用时 `orient_numbering` 返回 `None`，`number` 随即抛 `ValueError("numbering_failed")`；本层不产出半截命名结果。
- 指示氢、加氢前缀、位次省略标志都在本层定稿，Layer5 只做拼接与词形变化。

分派优先级（P-14 三层分派）：螺环与桥环的并列骨架候选裁决 → 保留模板固定编号（P-14.4(a)）→ 稠环外周编号（P-25.3.3）→ 通用 P-14.4 候选枚举。

输入 `parent` 的关键字段：`mol`、`chain`（L2 骨架原子顺序）、`scaffold_id` / `scaffold_match`（保留名模板）、`numbering_scaffold`（模板 locant 标签）、`principal_expression_facts`（主 FG 附着原子与类别，含 `multiplicity`）、`covered_principal_ids` / `principal_occurrences`、`bridged_node(s)` / `spiro_node(s)` / `fbs_node(s)`（并列骨架候选）、`double_bond(s)` / `triple_bond(s)`、`hydro_atoms`、`radical_c_idx`、`fused_tree`、`kind`、`n_carbons`；`substituents` 每项至少含 `attach_idx` / `en` / `kind`，稠合组分的稠合点带 `fusion` 标记（作编号取向依据，不当作取代基）。`component_numbering`（`fused_component_numbering` 内部合成 parent 时置真）令 `orient_numbering` 跳过指示氢收窄。

输出 numbered dict 的关键键：

- `parent`：`{**parent, "chain": 定向后的 chain}` 的浅拷贝，补写 `numbering_scaffold`（其 `labels` 为与链等长的位次标签，稠环为数字+字母位）、`hydro_prefix`（`(en, zh)` 加氢前缀，放弃时为空串对）、`hydro_fallback`（布尔，提示 L5 补写指示氢）、`indicated_h_locants`（指示氢位次字符串列表）、`indicated_h_forced`（布尔，保留名未隐含而须强制标注指示氢，或母体为生成式 HW 环）、`ind_h_carbon_ok`（布尔，母体是否为可写碳位指示氢的保留杂环名）、`lambda_locants`（`((位次, λ 键数|0, δ 双键数|0), …)`，需 λ/δ 标注的骨架原子）。
- `substituents`：每项附加 `locant`（无位次时取 0）。
- `fg_locants`：`[{kind, locants, omit}]` 稀疏列表，只产存在的项；`ene_locants` / `omit_ene_locant` / `yne_locants` / `omit_yne_locant` 给不饱和位次及省略标志。

数据流：`number(parent, substituents)` 先调 `orient_numbering` 得到定向后的 `chain`，用 `{**parent, "chain": chain}` 生成 `oriented`；`_with_locants` 给每个取代基补 `locant`；`_pack` 组装 `fg_locants` 与不饱和位次两组记录。随后在 `packed`（即 `oriented`）上依次处理加氢位、指示氢与加氢前缀，再补写 `lambda_locants` / `ind_h_carbon_ok`，最后返回 `result`。候选裁决在 `orient_numbering` 内写回调用方的 `parent`（`bridged_node` / `spiro_node` / `fbs_node` / `numbering_scaffold`），这些写入早于浅拷贝，故 `result["parent"]` 一并带上；`number` 自身写入的字段只落在该副本上。

实现约定：

- 失败即放弃：位次无法完整表达时不产出近似结果——`orient_numbering`、`fbs_numbering`、`bridged_numbering`、`number_fused_system` 返回 `None`，`hydro_prefix` 返回空串；`assign_cip` 与模板匹配用到的氢化骨架均经 `tools.memo` 按分子复用。
- 位次标签统一由 `numbering_scaffold.labels` 承载：稠环为数字+字母位（`4a` 一类），组分式螺环位次带撇号（`1'` 一类），其余为纯数字。

## 核心逻辑

### 定向编号引擎（分派）

`numbering_engine.orient_numbering(parent, substituents, *, float_hetero=False)` 是唯一的定向入口：`chain` 为空时返回 `None`，否则按固定优先级尝试下列路径，先命中者胜出。

1. `parent["fbs_nodes"]` 非 `None`（P-24.5~24.7 组分式螺环）：交给 `spiro_numbering.fbs_numbering`，无论成败均不回落。
2. `parent["spiro_nodes"]` 非 `None`（P-24 螺环）：交给 `spiro_numbering.spiro_numbering`，无论成败均不回落（`chain_fused` 判据对螺环恒真，会被 `_fused_numbering` 误吞）。
3. `parent["bridged_nodes"]` 非 `None`（P-23 桥环）：交给 `bridged_numbering.bridged_numbering`，无论成败均不回落（`chain_fused` 判据对桥环同样恒真）。
4. `_fixed_numbering`（P-14.4(a)）：母体命中保留模板且模板可全覆盖环系时按模板映射固定编号，直接返回而不进入候选枚举。
5. `_fused_numbering`（P-25.3.3）：多环稠合系统走几何优选取向 + 外周编号。
6. 通用候选枚举：单环与开链走 P-14.4 候选生成 + 逐层收窄；只带 `bridged_node`（无 `bridged_nodes`）时在此前返回 `None`（桥环 `chain` 非环序，旋转/翻转枚举无意义）。

`float_hetero=True` 时（由 `fused_component_numbering` 在 `shared` 非空时透传）跳过杂环指示氢 NH 位次的收窄，保留镜像候选，让 locant 1 由稠合原子定。

桥环与螺环编号共享的候选设施在 `numbering_engine`：`candidates(parent, key)` 读 `parent[f"{key}s"]` 并列候选表、缺失时回落为 `parent[key]` 的单元素表；`hetero_atoms` / `by_z` 取骨架杂原子并按原子序数分组；`alpha_locants(numbering, chain, substituents)` 给出 P-14.4(g) 的前修饰基位次元组；`pick_equivalent(nodes, feat)` 在并列候选特征全同才取首个，否则判为不可判定并返回 `None`；`node_feature_key(..., extra)` 是特征工厂，键为 `(descriptor, extra(node), 杂原子位次, 后缀位次, 取代基位次)`；`resolve_numbering(parent, substituents, node_key, feature_fn, ladder_fn)` 是桥环/螺环通用的裁决外壳：取候选 → 施加 `ladder_fn` 收窄阶梯 → `pick_equivalent` 定案 → 把选中节点写回 `parent[node_key]` → 按位次升序返回原子表（不可判定返回 `None`）。

### 通用候选收窄（P-14.4）

候选生成：`_ring_cands` 对环给出全部 `n` 个旋转与 `2n` 个翻转候选；开链给 `_numbered(chain)` 与 `_numbered(reversed(chain))` 两个方向。`_numbered` / `_to_chain` 在「原子顺序列表」与「`{原子: 位次}`」之间互转，收窄全程在后者上做。

收窄原语（L2、桥环与螺环编号亦复用）：

- `narrow(cands, key_fn, *, reverse=False, skip_none=False)`：保留 `key_fn` 键最小（`reverse=True` 取最大）的全部候选；候选不足 2 个时原样返回。`skip_none=True` 且存在 `None` 键时，表示该判据不适用，候选原样返回。
- `narrow_by_senior(cands, key_fn, heteros, by_z, *, skip_none=False)`：先按杂原子全集的位次集合收窄，再按 `P145_SENIOR` 逐元素收窄该元素原子的位次，实现 P-14.4(a)(b) 的元素序。

位次键：`_locant_set` 求某原子集在当前候选下的排序位次元组（原子集与候选无交集时返回 `None`）；`_edge_locants(pos, n, bonds)` 在已建的编号顺序映射上求一组多重键沿编号方向占据的边位置，对跨编号首尾的闭合边（seam）记 `n`，端点不沿编号相邻时返回 `None` 表示判据不适用；`_bond_locant_pairs(cand, bonds, doubles)` 只建一次编号顺序，同时给出「全部多重键」与「双键」两组边位次（P-14.4(e)）；`_fusion_side_index_cand(pos, atoms, n)` 求稠合键在编号序列中的侧序号（首起 1、收尾侧记 `n`，非键记 `n+1`），供 P-25.3.1.3 定母体组分取向。

收窄层级（通用路径，依次执行，候选剩 1 个即停）：

1. 杂环先行 `_narrow_hetero_ring`：位次 1 给 `P145_SENIOR` 引用序最先元素 → 杂原子集 → 元素序 → 指示氢 NH 位次最小化（P-22.2.2.1.3）；链上无杂原子时跳过。
2. 环系先做 `_hydro_indicated_atoms`：加氢位与指示氢位并集的位次最小化（P-14.4(b)）；链系此步后移至 (c) 之后。
3. `_principal_atoms`：principal 特征基团附着原子的位次集合最小化（P-14.4(c)）。
4. 链系补做 `_hydro_indicated_atoms`（P-14.4(b)(d)(e)(i)）。
5. 不饱和键位（P-14.4(e)）：`_bond_locant_pairs` 的键位次元组与双键位次元组一并最小化。
6. 稠合侧字母（P-25.3.1.3）：`_fusion_side_index_cand` 让母体组分的稠合侧落在尽可能靠前的位次，仅当有 `fusion` 标记的取代基时参与。
7. 取代基位次集合最小化（P-14.4(f)），`fusion` 标记者不计入。
8. 字母序平局：`_stem_loc_pairs` 产出 `(alpha_order_key(en), 位次)` 排序对，字母序最前的取代基得最低位次（P-14.4(g)）。
9. 立体平局（P-14.4(j)）：`_chain_rs_codes` 取母体链上的 CIP 代码，`_rs_locant_key` 把 `RS_HI`（R/M/r）位次排在前、`RS_LO`（S/P/s）位次排在后，小者优先。

`_hydro_indicated_atoms` 的构成：`parent["hydro_atoms"]` 非空时取饱和环位并入加氢位；无加氢前缀时，仅对可写碳位指示氢的保留杂环母体（经 `_ind_h_carbon_scaffold_ok`）取 `indicated_hydrogen_atoms` 的碳位，令其取最低位次（P-14.4(b)/P-58.2.1）；`component_numbering` 的稠合组分仿母体不参与。

特征提取辅助：`_is_ring(parent)` 按 `scaffold_id` 是否为真判断环系，决定候选来自环旋转还是链反向；`_unsat_bonds(parent)` 把标量字段 `double_bond` / `triple_bond` 与列表字段 `double_bonds` / `triple_bonds` 统一成「全部多重键、双键」两个端点对列表——P-14.4 需按「双键 / 全部不饱和键」分别比较。

`assign_cip(mol)` 经 `memo.by_mol("cip", ...)` 记忆化重算 CIP：隐式 H 手性碳先补显式 H，再 `AssignStereochemistry` + `rdCIPLabeler.AssignCIPLabels`，最后把 `_CIPCode` 回写原分子。`_chain_rs_codes` 只保留链上落在 `RS_HI` / `RS_LO` 的代码；`_rs_locant_key` 的位次取并行 `labels`（长度相符时），否则用链序号。

### 环系定向与固定编号

`_fixed_numbering`（P-14.4(a)）针对命中保留模板的母体：经 `_template_matches` 求模板 `_Q` 的全覆盖子结构匹配（`uniquify=False`），直接匹配为空时退到 `_Q_H` + `memo.by_mol("hydrogenated", _hydrogenated, mol)` 的氢化骨架匹配；每个匹配经 `standard_chain(sid, match)` 得标准链，只保留原子集与 `chain` 相同的映射。

仅有一个映射时直接返回；多个对称映射并存时用 `(稠合侧索引 P-25.3.1.3 / (b) 指示氢位次集 / (c) 后缀位次集 / (f) 前缀位次集 / (g) 引用序)` 五元键取最小，后缀集含 principal 附着原子与 `radical_c_idx`（自由价连接点与 principal 同属 (c) 类）。带 `fusion` 标记的取代基单列作取向依据、不计入前缀集；单环保留杂环（非并进更大稠合系者）另取 `saturated_ring_atoms` 作 (b) 指示氢位。三者皆空时取首个映射。位次键 `_locant_key_of` 优先用 `_STANDARD_LABELS` 标签（长度与标准链相符时，按 `standard_path` 标签而非链位置取值），否则用链序号。

`_component_labels(parent, chain)` 为稠合组分补并行 labels：优先取 `numbering_scaffold.labels`，其次 `_STANDARD_LABELS`，都缺失或长度不符时退化为纯数字。

### 稠环编号

`fused_orientation` 负责 P-25.3.2.3 的几何优选取向。

`horizontal_rows` 从每条共享边出发，沿 `opposite_bonds`（偶环 1 条对面边、奇环 2 条）经 `_extend_row` 递归扩展水平行，按长度降序返回全部行候选。`_layout` 按行摆放：首环用 `RING_TEMPLATES[n]`（正 n 边形，一条边竖直）落位，其余环以共享边对齐，`_opposite_side` 让新环落在前一环对面；行内奇环双侧融合时改用 `ring_shape_template` 的变形环模板，构造不出则该行弃用。行以外的环反复经 `_place_neighbor` 摆放（首选共享环对侧，与任一环重叠超 `OVERLAP_FRAC` 时改试对侧），摆不下则整行弃用（返回 `None`）。`_valid_deform_overlap` 用 `_ring_deform`（Kabsch 刚体+缩放拟合，循环移位与镜像取最优）检查每环相对模板的最大偏差不超过 `DEFORM_MAX`，并检查非共享环之间的重叠。

`preferred_orientations` 对每个候选行取 `flip` 镜像（`y → -y`），评分键 `(水平行环数, Q1 右上象限占比, -Q3 左下象限占比, 水平轴上方环数)` 取最大者，并列全部返回。象限占比由 `_quadrant_fractions` 以水平行中心为原点按环计数法累加；`_row_center` 对偶数行取中间一对环的共同键中点、奇数行取中心环质心。`Orientation` 是冻结 dataclass（`row` / `coords` / `quad`），`coord_dict()` 给出 `{原子: (x, y)}`；`preferred_orientation` 返回首个候选供前端使用。

`fused_numbering.number_fused_system(mol, rings, coords, sub_layers=None, alpha_subs=None)` 负责 P-25.3.3 外周编号。

- `_candidates` 枚举起点：`_top_rings` 取最上端水平行最右环（y 聚行后取 x 最大），起点环无起点原子时改用其稠合邻环；`_top_atoms` 取环内 y 最大的非稠合原子、优先与稠合原子相邻者，平局全保留；`_boundary_walk` 沿 `_exterior_edges`（只属一个环的边）行走外边界，非简单环返回空表，鞋带面积为正时整体反转以保证顺时针。`neighbors`（原子 → 环内邻居）与 `exterior`（只属一个环的边）与坐标无关，由 `number_fused_system` 按环集算一次传入。
- `_assign_labels` 生成标签：非稠合原子与稠合杂原子取下一数字，稠合碳原子取「紧邻前数字 + a/b/c 递增」的字母位，得到数字+字母混合 locant；`fused_atoms` 给出出现在 ≥2 个环的稠合原子。
- 基线收窄 (a)–(d)：杂原子集按 `P145_SENIOR` 逐元素最小化、稠合碳位次最小化、稠合杂原子位次最小化。
- `sub_layers` 分层收窄：哨兵 `INDICATED_H` 触发 `_as_indicated`（指示氢候选位由 `saturated_ring_atoms` 给出，候选数为奇数时才适用，且位次须全部落在链内）；其余层按位次集合最小化收窄镜像。
- 仍并列时按 `alpha_subs` 做 P-14.5 字母序收窄，最后按 `_rs_locant_key` 做 P-14.4(j) CIP 收窄。
- 返回值 `(chain, labels)`，调用方把 `numbering_scaffold` 记为 `{"scaffold_id": "fused", "labels": ...}`。

`_fused_numbering` 的准入与层序：`TRADITIONAL_NUMBERING_IDS` 中的骨架（anthracene、phenanthrene、acridine、carbazole、purine、xanthene、thioxanthene、cyclopenta[a]phenanthrene）按传统编号，从这条路径排除；`parent["bridged_node"]` 或 `parent["spiro_node"]` 非空时返回 `None`。走该路径的判据为「链原子全芳香、或命中已注册稠合模板（`get_spec(...).n_rings >= 2`）、或未注册稠环（SSSR 覆盖链的环数 ≥ 2）」，且取到的环系 SSSR 环数须 ≥ 2。`sub_layers` 层序为：游离价连接点（P-29）→ principal 特征基团（P-14.4(c)）→ `INDICATED_H`（P-25.3.3.1.2(f)，P-31.2.2 指示氢优先）→ 加氢位次 `hydro_atoms`（P-14.4(e)(i)）→ 取代基位次集（P-14.4(f)）→ 字母序平局（P-14.5）。

`fused_component_numbering(mol, scaffold_id, sub_rings, shared=None, sub_edges=None, *, side_letter=False)` 给单个稠合组分自身编号（P-25.3.3）。无 `scaffold_id` 的单环组分直接返回环原子表与纯数字标签；单环组分、或有 `_STANDARD_ORDERS` 固定编号的组分走 `orient_numbering`（合成 parent 置 `component_numbering=True`，`shared` 非空时置 `float_hetero=True`），标签由 `_component_labels` 补出；其余多环组分走 `preferred_orientations` + `number_fused_system`，把稠合点原子集作为单层 `sub_layers` 逐层最小化位次。`shared`（稠合点）以 `{"attach_idx": a, "fusion": side_letter}` 传入：`side_letter=True`（母体组分）走 P-25.3.1.3「稠合侧字母靠前」，附加组分仍最小化自身稠合位次。返回 `(chain, labels)`，失败返回 `(None, None)`。

`ring_geometry` 提供平面几何原语：`regular_polygon` 与 `RING_TEMPLATES`（3–98 元环）、`ring_shape_template`（5/7 元环变形模板）、`ring_cyclic`（以指定边为首的环序）、`rigid_fit` / `apply_rigid`（Kabsch 刚体拟合）、`clip_polygon`（Sutherland–Hodgman 裁剪）、`polygon_area` / `overlap_area` / `centroid`；阈值 `DEFORM_MAX` / `OVERLAP_FRAC` 为模块级常量，取向判定只依赖坐标与离散计数（象限值与环数取 0/0.25/0.5/1），无浮点累计尾差。

### 桥环裁决

`bridged_numbering.bridged_numbering(parent, substituents)` 只做裁决：von Baeyer 拆解与 P-23.2 拓扑收窄由 L2 完成，本层只在并列最优 `BridgedNode` 间按 P-23.3.2 与 P-14.4 定编号。

- 经 `resolve_numbering(parent, substituents, "bridged_node", lambda nd: nd.locant_pairs, _narrow_ladder)` 汇总：候选取自 `candidates(parent, "bridged_node")`，无候选或 `chain` 为空时返回 `None`。
- `_narrow_ladder` 依次施加判据（剩 1 个即停）：P-23.3.2.1 杂原子位次集合最低 → P-23.3.2.2 按元素序逐元素窄化（复用 `narrow_by_senior`，`by_z` 按原子序数分组）→ P-14.4(c) 后缀位次最低 → P-14.4(f) 取代基位次集合最低 → P-14.4(g) 字母序最前的取代基得最低位次（`alpha_locants`）→ P-31.1.4.2 残余平局（`_unsat_key`：复合位次数目最少 → 忽略括号后的较小端位次元组 → 全部位次元组）。P-23.3.2.2 的序列等价于 `P145_SENIOR` 去掉卤素——卤素一价，做不了骨架原子。
- 收窄后 `pick_equivalent` + `node_feature_key(..., lambda nd: nd.locant_pairs)`：剩余候选的渲染特征（`descriptor`、次级桥上标位次对、杂原子/后缀/取代基位次）全同才取首个，否则判为**不可判定**并返回 `None`。
- 命中后把选中节点写回 `parent["bridged_node"]`，返回其位次升序原子表供下游当 `chain` 用，两者须自洽。

### 螺旋环编号

`spiro_numbering.spiro_numbering(parent, substituents)` 处理 P-24 螺环：经 `resolve_numbering(parent, substituents, "spiro_node", lambda nd: nd.descriptor_superscripts, _narrow_ladder)` 汇总——`candidates(parent, "spiro_node")` 取候选、`_narrow_ladder` 收窄后由 `pick_equivalent` + `node_feature_key(..., lambda nd: nd.descriptor_superscripts)` 定案（并列候选描述符与螺原子上标不一致即不可判定），命中则把节点写回 `parent["spiro_node"]` 并返回位次升序原子表。

`_narrow_ladder` 依次施加（剩 1 个即停）：P-24.2.2.1 / P-24.2.3.1 螺原子位次集合最低 → P-24.2.2.2 / P-24.2.3.2 描述符数字取小（键为 `(descriptor, descriptor_superscripts)`）→ P-24.2.4.1.2 杂原子集合与逐元素窄化（`narrow_by_senior`）→ P-14.4(c) 主特征基团位次最低 → P-14.4(e) 双键位次最低（`_bond_locant_pairs`）→ P-14.4(f) 取代基位次集合最低 → P-14.4(g) 字母序最前的取代基位次最低。后缀判据先于不饱和判据，与 P-14.4 的 (c) 在 (e) 之前一致。

`spiro_numbering.fbs_numbering(parent, substituents)` 处理 P-24.5~24.7 组分式螺环：候选取自 `candidates(parent, "fbs_node")` 并取首个节点；任一组分无候选或组合不可判定时返回 `None`。

- `_shrink(comp, mol, parent, substituents)` 在单个组分内收窄编号候选：P-14.4(c) 主特征基团位次最低 → P-24.5.2 / P-24.5.4 螺稠合位次优先于 `a` 前缀位次 → 杂原子位次集合最低（`narrow_by_senior`，P-14.4(e) / P-22.2.2.1.3）→ P-14.4(f) 取代基位次集合最低。对称保留母体的自同构在此不可分辨，剩余并列交给 `_joint_pick`。
- `_joint_pick` 对组分候选做笛卡尔积，按 `_listing_key`（沿 `node.links` 展平、逐个螺稠合位次对取 `(a 组分位次, b 组分位次)` 并归一为 `locant_key`，P-24.6.1）取最小组合；位次缺失的组合跳过，并列取先出现者；组合总数超过 `_MAX_COMBOS`（4096）时返回 `None`。
- 拼装全局 `chain` 与 `labels`：第 k 个组分（`comp.index`）的位次带 k 个撇号（`APOSTROPHE`）；螺原子同属两组分，`seen` 守卫保证全局链每原子只留一次，取先列组分的位次。
- 命中后把 `node.components` 换成各组分选中的编号并写回 `parent["fbs_node"]`，同时写 `parent["numbering_scaffold"] = {"scaffold_id": "fused_bridged_spiro", "labels": ...}`；返回位次升序原子表。

`_fbs_locant_set(num, atoms)` 求指定原子集在组分候选下的位次集合（经 `locant_key` 归一），是组分内各层收窄共用的位次键。

### 位次计算

`locant_calc` 负责位次计算与打包。

- `locant_key(x)`：locant → `(数值, 字母尾+撇号)` 排序键，兼容数字位、字母位（`4a`）与撇号位（`1'`）；`locant_str_sort` 按此键排序。
- `atom_locant(chain, atom, facts)`（`_atom_locant` 的公开别名）：按 `facts["labels"]` 定位次（纯数字返回 int，字母位原样返回字符串），标签缺失或长度不符时退回链序号，`atom` 不在链内时返回 `None`；`_atom_locants` / `_typed_atom_locants` 把一组骨架原子映射为排序后的位次列表。
- `_typed_group_atoms(parent, group)` 取 `principal_expression_facts.attachment_atoms`，仅当其 `group_class.value` 等于 `group` 时非空——位次记录只对主官能团类别产生。
- `_occ_attachment` 把 occurrence 锚点映到骨架内附着原子（锚点全在骨架外时取骨架内邻居）；`_occurrence_locants` 逐 occurrence 求附着位次，只收 `covered_principal_ids` 中的 occurrence，同一原子承载多个同类 FG 时保留重数；`_expand_shared_locants` 对 `spec.locant_source == "attachment"` 的项，在 `facts.multiplicity` 大于已收位次数且逐 occurrence 位次个数恰等于 `multiplicity` 时，改用逐 occurrence 位次补回重复位次（偕二醇 `propane-2,2-diol`）。
- `_FG_GROUP` 记录 kind → `principal_expression_facts` 类别映射（键与 `FG_SPECS` 的 `fg` 同形：`alcohol` / `amine` / `ketone` / `thiol`）；`_omit_for` 据此把 FG 记录的省略判定接到 `omit_fg_locant`，`pos` 传 `True`，故实际判据落在「`carbocycle` 且非稠环」分支与 `n_carbons` 上，ketone 多原子时 `single=False`。
- `_FG_LOCANTS` 由 `FG_SPECS` 派生 `(spec.fg, spec)` 对，保证跨层 fg 一致性；`_fg_locants` 逐项取 `_typed_atom_locants(oriented, spec.fg)`，只产位次非空的记录 `{kind, locants, omit}`。
- 不饱和位次：`_as_bond_pairs` 把标量键字段与键列表字段统一成键对序列；`_bond_locant` 给单键位次，两端编号相邻取较小者，否则写复合位次 `x(y)`（P-31.1.4.2(1)）；`_bond_min_locs` / `_bond_locants` / `ene_locants` / `yne_locants` 取全部双键与三键的位次排序列表，`_unsat_locants` 打包 `ene_locants` / `omit_ene_locant` / `yne_locants` / `omit_yne_locant`。
- `omit_fg_locant(pos, n_carbons, parent, n_subs, *, single=True)`：`scaffold_id == "carbocycle"` 且非稠环且 `single` 时无取代则省略位次，否则位次为 1 且碳数不超过 2 时省略；`omit_unsat(n_carbons, kind, parent, *, triple=False)`：`kind == "heterane"` 的杂原子链（二核/三核的单一不饱和键）位次省略（diazene / triazene / disilyne，P-14.3.4.2(d)）；`kind == "alkane"` 的纯烃环无多双键时烯位次隐含省略，烯 ≤C2、炔 ≤C3 时位次 `1` 省略。
- `suffix_locant_set(numbered)` / `prefix_locant_set(numbered)`：P-44.1.1 的后缀位次集合与 P-45.2.2 的前缀位次集合（P-14.3.5），从 `numbered["parent"]` 的 `chain` 与 `numbering_scaffold.labels` 换算后按 `locant_key` 排序，无位次的项不入键；取 `_principal_atoms` 用函数内导入，避开与 `numbering_engine` 的循环依赖。

`tools.re.alpha_order_key(stem)` 是全库前缀引用顺序的统一排序键（P-14.5 字母数字序），键形为 `(字母序列, 前导位次元组)`：`alkyl_alpha_key` 循环剥除前导立体描述符组、方括号与位次集（含 `1H-` 指示氢前缀）得稳定词干，`_nonitalic_letters` 取该词干中的「非斜体罗马字母序列」——位次、连字符、括号、立体描述符一律不参与比较，词中斜体前缀 `tert-` / `sec-` 先用 `_MID_ITALIC_RE` 去掉（正则只匹配非字母边界后的 `tert`/`sec`，避免误伤词中的同形片段），斜体描述符 R/S/E/Z/H/N 为大写，`_LOWER_LETTER_RE` 只收小写字母，故同样不参与首轮比较；`_lead_locants` 取原词干首个非斜体罗马字母左端的位次集为次级键，无位次者最优先（P-14.5.4）。`numbering_engine` 的 `_stem_loc_pairs`、`_fixed_numbering`、`_fused_numbering`、`alpha_locants` 以及 Layer5 的各前缀排序都使用该键。

### 指示氢

`numbering_engine._nh_sites(heteros, mol)` 实现 P-14.4(b)/P-31.2.2，给出收窄杂环编号用的氮位集合：环内无 NH（`GetTotalNumHs() > 0`）时返回空表（位次不定）；否则取全部 NH 与全部无 H 的度为 3 的 N-取代氮，两者同为母体指示氢等价位，位次集并列，交由后续 (c)/(f) 裁决。由此避免把 `=N-` 型氮误当作母体氢化物的 NH 位。

`_narrow_hetero_ring(cands, mol, chain, float_hetero)` 把链上非碳原子收集为杂原子集，先令 `P145_SENIOR` 引用序最先的元素得位次 1，再走 `narrow_by_senior` 完成 P-22.2.2.1.3 的 (a)(b) 元素序收窄，`float_hetero=False` 时再按 `_nh_sites` 给出的 NH 位次集合收窄；杂环的元素序收窄先于主流程中的 principal 特征基团收窄。

`indicated_hydrogen` 模块实现 P-58.2.1 指示氢位次。

`saturated_ring_atoms(mol, ring_atoms, exclude=frozenset())` 返回「仅以单键连邻环原子且带氢」的饱和环位：入参先滤掉非环原子，无环位则返回空表；H 数取原始 `mol`（Kekulize 会给吡啶型 N 补隐式 H），键级取 `kekulized(mol)` 结果；`exclude` 中的位次由加氢前缀表达，不重复标指示氢；Kekulize 失败时芳香键保持 AROMATIC，该位静默放弃。`_flanked_by_exo_double` 护栏：环碳的两个环邻位都带环外杂原子多重键（如 1,3-二酮的 C2）时不列指示氢——不饱和度已由后缀确定（P-14.4）。

`indicated_hydrogen_atoms(mol, chain, exclude=frozenset(), extra=frozenset(), scaffold_id=None)` 是核心，合并 `saturated_ring_atoms` 结果与 `extra`（保留母体名未隐含、须显式标注的环位），按 `chain.index` 排序后输出环位原子表；`indicated_hydrogen(...)` 只是包装，把原子表按 `labels`（缺失用链序号）转成位次字符串列表。`is_hw(scaffold_id)` 判定母体是否为生成式 HW 杂单环（经 `layer2.hantzsch_widman.is_hw_scaffold`）。饱和位清空判据依次判定：

- 母体命中螺环保留名（`_is_spiro_scaffold`，经 `layer2.spiro_system.SPIRO_SCAFFOLDS`）时只留 `extra`——螺环名的 ene/yne 位次已完整表达氢化度。
- 链上原子构成单个环（`_is_monocycle`）且环内多重键数为 0（`_ring_double_bonds`）且无 `extra` 时清空——母体氢化物名已隐含全部 H。
- 单环且环内恰有一个多重键、母体未命中保留名模板（`_is_retained_scaffold`，经 `layer2.ring_scaffold.get_spec` 判定）且非 HW 名时清空——环己烯/环戊烯的氢位无歧义；HW 名按 mancude 词干读，须标指示氢。
- HW 名且无 `extra` 时，再由 `_forced_h_atom` 滤除价态强制的氢位（环邻位的有效价已容不下环内 π 键，如 1,3-二氧戊环的 C2；P-58.2.1）。

输出前还有一道互变异构冗余护栏：饱和位全为氮且多于一个时只保留最低位次。

`numbering.number(parent, substituents)` 是组装点：`orient_numbering` → `_with_locants` + `_pack` → 以 `packed` 的 `numbering_scaffold.labels` 与 `scaffold_id` 传入 `indicated_hydrogen`。加氢位 `hydro_atoms` 先经 `_fallback_hydro_atoms` 补未注册稠环推导（`fused_tree` 非空时取 `saturated_ring_atoms` 全部饱和环位，含芳香稠环上的 NH 位；碳子集原子数落在 `HYDRO_MULT_N` 时优先用碳子集），当缺省或 `hydro` 是推导结果的真子集且推导数落在 `HYDRO_MULT_N` 或其减一也落在时改用推导结果（奇数个饱和位亦可回退）。随后 `_lowest_extra_to_indicated`（最低位次落在加氢位时与指示氢位中位次最高者互换）与 `_odd_hydro_to_indicated`（加氢位个数为奇数时去掉最低者，P-51.1.1.4 / P-58.2.1.2）重分名次。被移出的位次若属模板未隐含的加氢位（`_extra_hydrogenated`，经 `layer2.ring_scaffold.extra_hydrogenated_atoms`）则置 `hydro_fallback`，提示 Layer5 补写指示氢，否则该段饱和度整段丢失。最后由 `_lambda_locants` 写 `lambda_locants`、由 `_ind_h_carbon_scaffold_ok` 写 `ind_h_carbon_ok`，并把 `indicated_h_locants` / `hydro_prefix` 一并落到 `packed`。

`_lambda_locants(mol, chain, labels)` 逐骨架原子产出 `(位次, λ 键数|0, δ 双键数|0)`：λ 只对 `HW_RING_LAMBDA_Z` 元素（硫族/氮氧的环内高价由 dioxo/oxide 前缀表达，名称里不复标 λ）且被标记者取 `bonding_number`；δ 计该原子直接相连的骨架内连续双键数，仅 ≥2 时才记（P-25.7.2）。`_ind_h_carbon_scaffold_ok(scaffold_id)` 判定保留母体名是否可写碳位指示氢（`naming_class` 不为 `benzodioxole` / `adamantane` / `mono_carbo`——这几种的饱和位由保留名本身表达）。

`hydro_prefix(chain, labels, hydro_atoms)` 产出加氢前缀 `(en, zh)`：链或加氢位为空、数量不落在 `HYDRO_MULT_N` 时返回空串；完全氢化（链上原子全在加氢位内）时省略全部位次（P-14.3.4.5），否则列出位次；加氢原子不全在编号链内时返回空串，`number` 随即丢弃加氢位表达并整体退回指示氢，不产出半截名。`_extra_indicated` 经 `layer2.ring_scaffold.extra_indicated_atoms` 取保留模板同位无 H 而分子有 H 的芳香杂环原子；`indicated_h_forced` 取 `bool(extra)` 或母体为生成式 HW 环。`exclude`（`hydro`）与 `extra` 是方向相反的两类修正：前者把已由加氢前缀表达的位次排除，后者把保留母体名未隐含的位（如稠合杂芳环的 NH）补入。

## 与其他层的契约

对外入口：`numbering.number` / `hydro_prefix`；`numbering_engine.orient_numbering` / `fused_component_numbering` / `resolve_numbering` / `narrow` / `narrow_by_senior` / `assign_cip` / `_principal_atoms` / `_ring_cands` / `_narrow_hetero_ring` / `candidates` / `hetero_atoms` / `by_z` / `alpha_locants` / `pick_equivalent` / `node_feature_key`；`bridged_numbering.bridged_numbering`；`spiro_numbering.spiro_numbering` / `fbs_numbering`；`fused_numbering.number_fused_system` / `fused_atoms`；`fused_orientation.preferred_orientations` / `preferred_orientation`；`indicated_hydrogen.indicated_hydrogen` / `indicated_hydrogen_atoms` / `saturated_ring_atoms` / `is_hw`；`locant_calc.locant_key` / `locant_str_sort` / `atom_locant` / `omit_fg_locant` / `omit_unsat` / `suffix_locant_set` / `prefix_locant_set`；`tools.re.alpha_order_key`。

L2 侧：

- 入口 `numbering.number` 由 `namer._assemble_candidate` 调用，取代基先经 `namer._subs_for_numbering` 筛选与重挂；`ValueError("numbering_failed")` 由该调用点捕获，候选判为失败。`namer` 还以 `suffix_locant_set` / `prefix_locant_set` 写入候选 `meta` 的 `p44_1_1_key` / `p45_2_2_key`（P-44.1.1 / P-45.2.2）作比较键。
- `layer2.fused_system` 复用 `narrow` 做 (a)–(f) 逐准则收窄，并调用 `fused_atoms` / `number_fused_system` / `preferred_orientations` / `locant_key` 完成 (g)–(j) 依赖编号的收窄。
- `layer2.bridged_system` 与 `layer2.parent_skeleton` 复用 `narrow` 收窄并列候选；`layer2.spiro_system` 复用 `_ring_cands` 与 `_narrow_hetero_ring` 生成单环组分编号候选，并经 `fused_component_numbering` 委托多环组分编号。
- `layer2.ring_scaffold` 提供 `get_spec`、`_Q` / `_Q_H`、`_hydrogenated`、`standard_chain`、`_STANDARD_LABELS` / `_STANDARD_ORDERS`、`numbering_scaffold_facts`、`extra_indicated_atoms`、`extra_hydrogenated_atoms`，是本层固定编号、模板匹配与指示氢护栏的数据来源。
- `layer2.hantzsch_widman` 提供 `is_hw_scaffold`（生成式 HW 杂单环判定，决定指示氢清空与强制标注）、`HW_RING_LAMBDA_Z` / `effective_valence`（λ 标注与价态强制氢判定）；`tools.lambda_notation` 的 `bonding_number` / `is_lambda_marked` 供 `_lambda_locants` 取 λ 键数。
- `layer2.principal_expression` 写入 `principal_expression_facts`、`hydro_atoms` 与 `bridged_node(s)` / `spiro_node(s)` / `fbs_node(s)`，供本层 P-14.4(c)、P-31.2.2 与骨架候选裁决使用（见 [[architecture/layer2-parent-selector]]）。

L5 侧：

- `layer5.assembler` 读 `parent.indicated_h_locants` / `indicated_h_forced` / `ind_h_carbon_ok` / `hydro_prefix` / `hydro_fallback` / `numbering_scaffold.labels`，由 `join_hydro_prefix` 决定加氢前缀与指示氢的拼接形态（`ind_h_carbon_ok` 作保留杂环母体碳位指示氢的补写门）；`_lambda_prefix` 读 `lambda_locants` 渲染 λ/δ 前缀（P-25.6 / P-25.7.2）。
- `layer5.assembler._ensure_parent_stem` 按 `fbs_node` → `spiro_node` → `bridged_node` → 稠环 → 生成式词干的次序注入母体词干，描述符与上标位次取自本层选中的同一节点（见 [[architecture/layer5-name-assembly]]）。
- `layer5.chain_engine` 按 `fg_locants` 的 `kind` 查找主官能团位次，读 `omit_ene_locant` / `omit_yne_locant` / `ene_locants` / `yne_locants` 决定段式后缀。
- `layer5.assembler_prefixes` 用 `locant_str_sort` 与 `alpha_order_key` 排序前缀引用顺序；`layer5.stereo` 用 `locant_key` 与 `assign_cip`；`layer5.fused_namer` 用 `fused_component_numbering`；`layer5.spiro_namer` 用 `locant_key` 与 `spiro_node` / `fbs_node` 生成螺环母体名。

字段口径见 [[reference/core-data-contracts]]。
