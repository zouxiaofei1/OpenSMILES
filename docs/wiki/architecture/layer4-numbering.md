# Layer4: 编号与位次 (Numbering & Locants)

源文件 10 个 `.py`（约 1564 行）｜对外接口 `numbering.number`

## 概述

Layer4 承接 Layer2 选定的母体 `parent` 与 Layer3 提取的取代基 `substituents`，按 P-14.4 / P-25.3.3 / P-22.2.2 给母体骨架定向编号，并把编号结果物化为 Layer5 可直接消费的位次数据。

职责边界：

- 只决定骨架的「原子顺序 / 位次」，不生成名称文本、不选母体、不改动取代基集合。
- 编号不适用时 `orient_numbering` 返回 `None`，由上层回退，本层不产出半截命名结果。
- 指示氢、加氢前缀、位次省略标志都在本层定稿，Layer5 只做拼接与词形变化。

输入 `parent` 的关键字段：`mol`（RDKit 分子）、`chain`（L2 给出的母体骨架原子顺序）、`scaffold_id` 与 `scaffold_match`（保留名母体模板）、`numbering_scaffold`（模板 locant 标签）、`principal_expression_facts`（P-14.4(c) 主特征基团附着原子）、`double_bond` / `double_bonds` / `triple_bond` / `triple_bonds`（不饱和键）、`hydro_atoms`（L2 判定的加氢环位）、`radical_c_idx`（游离价连接点）、`fused_tree`、`ring_attach_idx`、`kind`、`n_carbons`。`substituents` 每项至少含 `attach_idx`（附着原子）、`en`（英文前缀）、`kind`。

输出 numbered dict 的关键键：

- `parent`：定向后的 parent，`chain` 已换成定向结果，并就地补写 `indicated_h_locants`、`indicated_h_forced`、`hydro_prefix`、`numbering_scaffold`。
- `substituents`：每项附加 `locant`（无位次时取 0）。
- `fg_locants`：`[{kind, locants, omit}]` 稀疏列表，只产存在的项。
- `ene_locants` / `omit_ene_locant` / `yne_locants` / `omit_yne_locant`：不饱和位次及省略标志。
- `parent.numbering_scaffold.labels`：与链等长的位次标签（稠环为数字+字母位）。

数据流：`number(parent, substituents)` 先调 `orient_numbering` 得到定向后的 `chain`，用 `{**parent, "chain": chain}` 生成 `oriented`；`locant_calc._with_locants` 给每个取代基补 `locant`；`locant_calc._pack` 组装 `fg_locants` 与不饱和位次两组记录。随后在 `packed`（即 `oriented`）上依次处理加氢位与指示氢，最后返回 `result`。`result["parent"]` 与传入的 `parent` 是同一 dict 对象，本层对 `parent` 的字段写入对调用方可见。

文件分工：`numbering.py`（入口与加氢前缀）、`numbering_engine.py`（候选枚举、收窄原语、三层分派、稠合组分编号）、`fused_numbering.py`（P-25.3.3 外周编号）、`fused_orientation.py`（优选取向）、`ring_geometry.py`（平面几何原语）、`indicated_hydrogen.py`（指示氢位次）、`locant_calc.py`（位次计算与打包）、`omit_locants.py`（位次省略规则）、`candidate_keys.py`（并列候选比较键）、`__init__.py`。

实现约定：

- 记忆化：`numbering_engine.assign_cip` 经 `tools.memo.by_mol` 按分子复用确定性中间结果，命名期内同一分子不重复重算 CIP。
- 失败即放弃：位次无法完整表达时不产出近似结果——`orient_numbering` 与 `number_fused_system` 返回 `None`，`hydro_prefix` 返回空串，`saturated_ring_atoms` 在 Kekulize 失败时跳过该位。
- 就地写入：`number` 不复制 `parent`，`result["parent"]` 与传入的 `parent` 为同一对象，本层写入的字段对调用方可见。
- 位次标签统一由 `numbering_scaffold.labels` 承载：稠环为数字+字母位（`4a` 一类），其余为纯数字；凡涉及位次比较与排序的位置一律经 `locant_key` 归一。

## 三层分派

`numbering_engine.orient_numbering(parent, substituents, *, float_hetero=False)` 是唯一的定向入口，按固定优先级尝试三条路径，先命中者胜出。

1. `_fixed_numbering`（P-14.4(a)）：母体命中保留模板且模板可全覆盖环系时，按模板的 `standard_chain` 映射固定编号，直接返回而不进入候选枚举。多个对称映射并存时，用 `(c) 后缀位次集 / (f) 前缀位次集 / (g) 引用序` 三元键取最小；无后缀也无前缀时直接取首个映射。
2. `_fused_numbering`（P-25.3.3）：多环稠合系统走几何优选取向 + 外周编号。`TRADITIONAL_NUMBERING_IDS` 中的骨架（anthracene、phenanthrene、acridine、carbazole、purine、xanthene、thioxanthene、cyclopenta[a]phenanthrene）按传统编号，从这条路径排除。
3. 通用候选枚举：单环与开链走 P-14.4 候选生成 + 逐层收窄。

`float_hetero=True` 时（由 `fused_component_numbering` 在对称杂环上透传）跳过杂环指示氢 NH 位次的收窄，保留镜像候选，让 locant 1 由稠合原子定。

## 通用编号引擎

`numbering_engine` 负责候选枚举与 P-14.4 逐层收窄。

候选生成：`_ring_cands` 对环给出全部 `n` 个旋转与 `2n` 个翻转候选；开链给 `_numbered(chain)` 与 `_numbered(reversed(chain))` 两个方向。`_numbered` / `_to_chain` 在「原子顺序列表」与「`{原子: 位次}`」之间互转，收窄全程在后者上做。

收窄原语（公共原语，L2 亦复用）：

- `narrow(cands, key_fn, *, reverse=False, skip_none=False)`：保留 `key_fn` 键最小（`reverse=True` 取最大）的全部候选。`skip_none=True` 且存在 `None` 键时，表示该判据不适用，候选原样返回。
- `narrow_by_senior(cands, key_fn, heteros, by_z, *, skip_none=False)`：先按杂原子全集的位次集合收窄，再按 `P145_SENIOR` 逐元素收窄该元素原子的位次，实现 P-14.4(a)(b) 的元素序。

位次键：`_locant_set` 求某原子集在当前候选下的排序位次元组；`_bond_locants` 求多重键沿编号方向占据的边位置，对跨编号首尾的闭合边（seam）记 `n`，端点不沿编号相邻时返回 `None` 表示判据不适用。

收窄层级（通用路径，依次执行，候选剩 1 个即停）：

1. 杂环先行 `_narrow_hetero_ring`：杂原子集 → 元素序 → 指示氢 NH 位次最小化（P-22.2.2.1.3）。
2. `_principal_atoms`：principal 特征基团附着原子的位次集合最小化（P-14.4(c)）。
3. 不饱和键位（P-14.4(e)）：键位次元组与双键位次元组依次最小化。
4. 取代基位次集合最小化（P-14.4(f)）。
5. 字母序平局：`_stem_loc_pairs` 产出 `(alpha_order_key(en), 位次)` 排序对，字母序最前的取代基得最低位次（P-14.5）。
6. 立体平局（P-14.4(j)）：`_chain_rs_codes` 取母体链上的 CIP 代码，`_rs_locant_key` 把 `RS_HI`（R/M/r）位次排在前、`RS_LO`（S/P/s）位次排在后，小者优先。

特征提取辅助：`_is_ring(parent)` 按 `scaffold_id` 是否为真判断环系，决定候选来自环旋转还是链反向；`_unsat_bonds(parent)` 把标量字段 `double_bond` / `triple_bond` 与列表字段 `double_bonds` / `triple_bonds` 统一成「全部多重键、双键」两个端点对列表，单数键亦进「全部」列表——P-14.4 自身需按「双键 / 全部不饱和键」分别比较。

CIP 由 `assign_cip` 经 `memo.by_mol("cip", ...)` 记忆化重算：隐式 H 手性碳先补显式 H，再 `AssignStereochemistry` + `rdCIPLabeler.AssignCIPLabels`，最后只把 `_CIPCode` 回写原分子。`_chain_rs_codes` 只保留链上落在 `RS_HI` / `RS_LO` 的代码；`_rs_locant_key` 的位次取值优先用并行 `labels`（长度与链相符时），否则用链序号。

`_component_labels` 为稠合组分补并行 labels：优先取 `numbering_scaffold.labels`，其次 `_STANDARD_LABELS`，都缺失时退化为纯数字。

## 稠合系统编号

`fused_orientation` 负责 P-25.3.2.3 的几何优选取向。

`horizontal_rows` 从每条共享边出发，沿 `opposite_bonds`（偶环 1 条对面边、奇环 2 条）向两侧递归扩展水平行，按长度降序返回全部行候选。`_layout` 按行摆放：首环用 `RING_TEMPLATES[n]`（正 n 边形，一条边竖直）落位，其余环以共享边对齐，`_opposite_side` 让新环落在前一环的对面；行内奇环双侧融合时改用 `ring_shape_template` 的变形环模板，构造不出则该行弃用。行以外的环由 `_place_neighbor` 递归摆放，重叠面积超阈值时改试对侧。`_valid_deform_overlap` 用 `rigid_fit`（Kabsch 2D 刚体+缩放拟合）检查每环相对模板的最大偏差不超过 `DEFORM_MAX`，并检查非共享环之间重叠不超过 `OVERLAP_FRAC`。

`preferred_orientations` 对每个候选行取 `flip` 镜像，评分键为 `(水平行环数, Q1 右上象限占比, -Q3 左下象限占比, 水平轴上方环数)`，取最大者，并列全部返回。象限占比由 `_quadrant_fractions` 以水平行中心为原点、按环计数法累加；`_row_center` 对偶数行取中间一对环的共同键中点、奇数行取中心环质心。`Orientation` 是冻结 dataclass，`coord_dict()` 给出 `{原子: (x, y)}`。

`fused_numbering.number_fused_system(mol, rings, coords, sub_layers=None, alpha_subs=None)` 负责 P-25.3.3 外周编号。

- `_candidates` 枚举起点：`_top_rings` 取最上端水平行最右环（y 聚行后取 x 最大），`_top_atoms` 在环内取 y 最大的非稠合原子、优先与稠合原子相邻者，平局全保留；`_boundary_walk` 沿 `_exterior_edges` 单闭环行走外边界，鞋带面积为正时整体反转以保证顺时针。
- `_assign_labels` 生成标签：非稠合原子与稠合杂原子取下一数字，稠合碳原子取「紧邻前数字 + a/b/c 递增」的字母位，得到数字+字母混合 locant。
- 基线收窄 (a)–(d)：杂原子集按 `P145_SENIOR` 逐元素最小化、稠合碳位次最小化、稠合杂原子位次最小化。
- `sub_layers` 分层收窄：哨兵 `INDICATED_H` 触发 `_as_indicated`（指示氢候选位由 `saturated_ring_atoms` 给出，候选数为奇数时才适用）；其余层按位次集合最小化收窄镜像。
- 仍并列时按 `alpha_subs` 做 P-14.5 字母序收窄，最后按 `_rs_locant_key` 做 P-14.4(j) CIP 收窄。
- 返回值 `(chain, labels)`，调用方把 `numbering_scaffold` 记为 `{"scaffold_id": "fused", "labels": ...}`。

`_fused_numbering` 构造 `sub_layers` 的层序为：游离价连接点（P-29）→ principal 特征基团（P-14.4(c)）→ `INDICATED_H`（P-25.3.3.1.2(f)，P-31.2.2 指示氢优先）→ 加氢位次 `hydro_atoms`（P-14.4(e)(i)，先于可分离前缀）→ 取代基位次集（P-14.4(f)）→ 字母序平局（P-14.5）。是否走该路径由「链原子全芳香、或命中已注册稠合模板（`get_spec(...).n_rings >= 2`）、或未注册稠环（SSSR 覆盖链的环数 ≥ 2）」三者之一决定。

`fused_component_numbering(mol, scaffold_id, sub_rings, shared=None, sub_edges=None)` 给单个稠合组分自身编号（P-25.4 / P-25.3.3），稠合点当作取代基处理。单环组分、或有 `_STANDARD_ORDERS` 固定编号的组分走 `orient_numbering`（`shared` 非空时置 `float_hetero=True`，让对称杂环保留镜像、由稠合原子定 locant 1），标签由 `_component_labels` 补出；其余多环组分走 `preferred_orientations` + `number_fused_system`，把稠合点原子集作为单层 `sub_layers` 逐层最小化位次（P-25.3.1.3）。返回 `(chain, labels)`，失败返回 `(None, None)`。

`ring_geometry` 提供平面几何原语：`regular_polygon` 与 `RING_TEMPLATES`（3–19 元环）、`ring_shape_template`（5/7 元环变形模板）、`ring_cyclic`（求以指定边为首的环序）、`rigid_fit` / `apply_rigid`（Kabsch 刚体拟合）、`clip_polygon`（Sutherland–Hodgman 裁剪）、`polygon_area` / `overlap_area` / `centroid`。模板与阈值 `DEFORM_MAX` / `OVERLAP_FRAC` 是本模块的模块级常量，取向判定完全基于坐标与离散计数，评分键的象限值与环数为 0/0.25/0.5/1 的离散值，不引入浮点累计尾差。

## 杂环编号与指示氢

`numbering_engine._nh_sites(heteros, mol)` 实现 P-22.2.2.1.4，收窄杂环编号用的氮位集合：既取环内可带 H 的氮（`GetTotalNumHs() > 0`），也取无 H 但度为 3 的 N-取代氮，后者需满足「环内已有带 H 的氮」或「该氮在五元环内（`IsInRingSize(5)`）」之一。由此避免把 `=N-` 型氮误当作母体氢化物的 NH 位。

`indicated_hydrogen` 模块实现 P-58.2.1 指示氢位次。

`saturated_ring_atoms(mol, ring_atoms, exclude=frozenset())` 返回「仅以单键连邻环原子且带氢」的饱和环位：H 数取原始 `mol`（Kekulize 会给吡啶型 N 补隐式 H），键级取 `kekulized(mol)` 结果；`exclude` 中的位次由加氢前缀表达，该位不重复标指示氢；Kekulize 失败时芳香键保持 AROMATIC，该位静默放弃。

`_flanked_by_exo_double` 护栏：环碳的两个环邻位都带环外杂原子多重键（如 1,3-二酮的 C2）时不列指示氢——不饱和度已由后缀确定（P-14.4）。

`_is_monocycle` + `_ring_double_bonds == 1` + `_is_retained_scaffold(scaffold_id)` 护栏：链上原子构成单环、环内恰有一个多重键、且母体未命中保留名模板（`_is_retained_scaffold` 经 `layer2.ring_scaffold.get_spec` 判定）时，饱和位集清空——环己烯/环戊烯的氢位无歧义。

互变异构冗余护栏：饱和位全为氮且多于一个时只保留最低位次。

`_narrow_hetero_ring(cands, mol, chain, float_hetero)` 是杂环的分派差异点：把链上非碳原子收集为杂原子集并按原子序分组，先走 `narrow_by_senior` 完成 P-22.2.2.1.3 的 (a)(b) 元素序收窄；`float_hetero=False` 时再按 `_nh_sites` 给出的 NH 位次集合收窄。杂环的元素序收窄先于 `numbering_engine` 主流程中的 principal 特征基团收窄。

`indicated_hydrogen(mol, chain, labels=None, exclude=frozenset(), extra=frozenset(), scaffold_id=None)` 合并 `saturated_ring_atoms` 结果与 `extra`（保留母体名未隐含、须显式标注的环位，位次由 `chain.index` 排序），输出位次字符串列表；位次取 `labels`，缺失时用链序号。`exclude` 与 `extra` 是方向相反的两类修正：前者把已由加氢前缀表达的位次排除，后者把保留母体名未隐含的位（如稠合杂芳环的 NH）补入。`number` 中的 `indicated_h_forced` 即 `bool(extra)`，Layer5 据此判断是否必须强制写出指示氢而不依赖前缀省略规则。

`numbering.number(parent, substituents)` 是组装点，流程为：`orient_numbering` → `_pack` → 取 `packed` 的 `scaffold_id` 作 `scaffold_id` 传入 `indicated_hydrogen`。加氢位 `hydro_atoms` 先经 `_fallback_hydro_atoms` 补未注册稠环推导（取环内单键连邻环的饱和 sp3 位，碳原子数落在 `HYDRO_MULT_N` 时优先用碳子集），再经 `_lowest_extra_to_indicated`（最低加氢位次与指示氢位中位次最高者互换）与 `_odd_hydro_to_indicated`（加氢位为奇数个时最低者改用指示氢，P-51.1.1.4 / P-58.2.1.2）重分名次。`hydro_prefix(chain, labels, hydro_atoms)` 产出加氢前缀 `(en, zh)`：完全氢化时省略全部位次（P-14.3.4.5），否则列出位次；数量不落在 `HYDRO_MULT_N`、或加氢原子不全在编号链内时返回空串，`number` 随即整体退回指示氢表达而非产出半截名。`_extra_indicated` 经 `layer2.ring_scaffold.extra_indicated_atoms` 取保留模板同位无 H 而分子有 H 的芳香杂环原子，写入 `indicated_h_forced`。

## 位次省略与字母数字序

`locant_calc` 负责位次计算与省略判据。

- `atom_locant(chain, atom, facts)`（`_atom_locant` 的公开别名）：按 `facts["labels"]` 定位次（数字位返回 int，字母位原样返回字符串），标签缺失或长度不符时退回链序号。
- `locant_key(x)`：locant → `(数值, 字母尾)` 排序键，兼容数字与字母位混合；`locant_str_sort` 按此键排序 locant 集合。
- `_FG_GROUP`：记录 kind → `principal_expression_facts` 类别映射，键与 `FG_SPECS` 的 `fg` 同形：`{"alcohol": "alcohol", "amine": "amine", "ketone": "ketone", "thiol": "thiol"}`。`_omit_for` 据此把 FG 记录的省略判定接到 `omit_locants.omit_fg_locant`（ketone 多原子时 `single=False`）。`_FG_LOCANTS` 由 `FG_SPECS` 派生 `(记录 kind, spec)` 对，保证跨层的 fg 一致性。
- `_locants_for` 按 `spec.locant_source` 取位次：`anchor_field` 走 `_anchor_field_locants`（如 `radical_c_idx`），`attachment_exocyclic` 仅在该 FG 确以环外方式表达（`_exocyclic_only`）时取附着原子，否则走 `_typed_atom_locants`。
- `_unsat_locants` / `ene_locants` / `yne_locants` / `_bond_locants` / `_bond_min_locs`：取每根不饱和键较小端点的位次，排序后打包，并附 `omit_ene_locant` / `omit_yne_locant`。
- `omit_locants.omit_fg_locant(pos, n_carbons, parent, n_subs, *, single=True)`：`carbocycle` 且非稠环的环单 FG 无取代时省略位次；否则位次为 1 且碳数不超过 2 时省略。
- `omit_locants.omit_unsat(n_carbons, kind, parent, *, triple=False)`：纯烃环无多双键时烯位次隐含省略；烯 ≤C2、炔 ≤C3 时位次 `1` 省略（P-14.3.4.2(d)）。
- `_pack(oriented, substituents)` 组装最终 dict：`parent` / `substituents` / `fg_locants` / 不饱和位次四组键。

`tools.re.alpha_order_key(stem)` 是全库前缀引用顺序的统一排序键（P-14.5 字母数字序）。实现分两步：先 `_strip_n_prefix` 剥掉 `N-` / `N,` 前缀，再调 `alkyl_alpha_key` 循环剥除斜体前缀、括号、位次集（含 `1H-` 指示氢前缀）与立体描述符组，得到稳定词干；随后

1. `_nonitalic_letters` 取该词干中的「非斜体罗马字母序列」——位次、连字符、括号、立体描述符一律不参与比较，词中斜体前缀 `tert-` / `sec-` 先用 `_MID_ITALIC_RE` 去掉（正则只匹配非字母边界后的 `tert`/`sec`，避免误伤词中的同形片段）；斜体描述符 R/S/E/Z/H/N 为大写，`_LOWER_LETTER_RE` 只收小写字母，故同样不参与首轮比较。
2. `_lead_locants` 取首个非斜体罗马字母左端的位次集为次级键，无位次者最优先（P-14.5.4 平局判据）。

键形为 `(字母序列, 前导位次元组)`。`numbering_engine` 的 `_stem_loc_pairs`、`_fixed_numbering`、`_fused_numbering` 以及 Layer5 的各前缀排序都使用该键。

## 跨层可复用入口

```
numbering.number(parent, substituents) -> dict
numbering.hydro_prefix(chain, labels, hydro_atoms) -> (str, str)
numbering_engine.orient_numbering(parent, substituents, *, float_hetero=False) -> list[int] | None
numbering_engine.fused_component_numbering(mol, scaffold_id, sub_rings, shared=None, sub_edges=None)
numbering_engine.narrow(cands, key_fn, *, reverse=False, skip_none=False) -> list
numbering_engine.narrow_by_senior(cands, key_fn, heteros, by_z, *, skip_none=False) -> list
numbering_engine.assign_cip(mol) -> None
fused_numbering.number_fused_system(mol, rings, coords, sub_layers=None, alpha_subs=None)
fused_numbering.fused_atoms(rings) -> set[int]
fused_orientation.preferred_orientations(mol, rings, fusion_edges) -> list[Orientation]
indicated_hydrogen.indicated_hydrogen(mol, chain, labels=None, exclude=..., extra=..., scaffold_id=None)
indicated_hydrogen.saturated_ring_atoms(mol, ring_atoms, exclude=frozenset()) -> list[int]
locant_calc.locant_key(x) / locant_str_sort(locs) / atom_locant(chain, atom, facts) / _pack(...)
candidate_keys.suffix_locant_set(numbered) / prefix_locant_set(numbered) -> tuple
omit_locants.omit_fg_locant(...) / omit_unsat(...) -> bool
tools.re.alpha_order_key(stem) -> tuple
```

## 与 L2/L5 的耦合

与本层的接口面：

- 入口 `numbering.number(parent, substituents)` 由 `namer._assemble_candidate` 调用，取代基先经 `namer._subs_for_numbering` 筛选与重挂。
- `numbering_engine.orient_numbering`、`narrow`、`narrow_by_senior`、`assign_cip`、`_principal_atoms`、`_chain_rs_codes`、`_rs_locant_key`、`fused_component_numbering` 均为跨层复用的公共原语。
- `fused_numbering.fused_atoms`、`number_fused_system`、`fused_orientation.preferred_orientations` 供 L2 与 L5 直接调用。
- `candidate_keys.suffix_locant_set` / `prefix_locant_set` 给出 P-44.1.1 / P-45.2.2 的并列候选比较键，`namer` 把它写入结果的 `meta`。

L2 侧：

- `layer2.fused_system` 复用 `narrow` 做 (a)–(f) 逐准则收窄，并调用 `fused_atoms` / `number_fused_system` / `preferred_orientations` / `locant_key` 完成 (g)–(j) 依赖编号的收窄。
- `layer2.ring_scaffold` 提供 `get_spec`、`_Q` / `_Q_H`、`standard_chain`、`_STANDARD_LABELS` / `_STANDARD_ORDERS`、`numbering_scaffold_facts`、`extra_indicated_atoms`，是本层固定编号与指示氢护栏的数据来源。
- `layer2.kind_registry` 经 `numbering_scaffold_facts` 给 parent 挂 `numbering_scaffold`。
- `layer2.principal_expression` 写入 `principal_expression_facts` 与 `hydro_atoms`，供本层 P-14.4(c) 与 P-31.2.2 使用。

L5 侧：

- `layer5.assembler` 读 `parent.indicated_h_locants` / `indicated_h_forced` / `hydro_prefix` / `numbering_scaffold.labels`，由 `join_hydro_prefix` 决定加氢前缀与指示氢的拼接形态。
- `layer5.chain_engine` 按 `fg_locants` 的 `kind` 查找主官能团位次，读 `omit_ene_locant` / `omit_yne_locant` / `ene_locants` / `yne_locants` 决定段式后缀。
- `layer5.assembler_prefixes` 用 `locant_str_sort` 与 `alpha_order_key` 排序前缀引用顺序；`layer5.stereo` 用 `locant_key` 与 `assign_cip`；`layer5.fused_namer` 用 `fused_component_numbering` 给稠合组分编号。
