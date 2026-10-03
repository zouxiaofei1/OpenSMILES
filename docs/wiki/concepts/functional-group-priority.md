# 官能团优先级体系 (Functional Group Priority)

本文说明 NamePredict 中官能团优先级的声明处、跨层流转路径与仲裁规则。优先级数据以 `src/namepredict/layer1/fg_registry.py` 为唯一事实来源，L2 与 L5 只消费、不重新登记。

## 优先级从哪来

`fg_registry.FG_SPECS` 是全部官能团类别元数据的唯一声明处，每个条目是一个冻结 dataclass `FgSpec`。它的字段含义如下。

| 字段 | 含义 |
| --- | --- |
| `fg` | 官能团类别字符串（`FunctionalGroupClass` 值），同时用作 L1 analyzer 的列表 key |
| `p41` | P-41 表 4.1 的主官能团类号，`0` 表示不作主官能团 |
| `path` | P-43 优先级路径，用于同类号内部再排序 |
| `expr` | 表达类型（`suffix` / `prefix_only` / `legacy_compat`） |
| `anchors` | occurrence payload 中的锚点 key，空表示不收集锚点 |
| `parent_anchor_fields` | 母体锚点字段名（单数, 复数） |
| `locant_source` | 位次原子来源：`attachment` / `attachment_exocyclic` / `anchor_field` |

`FG_SPECS` 的全部条目及其 `p41` 如下。

| 顺序 | `fg` | `p41` | `path` | 备注 |
| --- | --- | --- | --- | --- |
| 1 | `radical` | 1 | `()` | 自由基，与酰基同级；`parent_anchor_fields=("radical_c_idx",)` |
| 2 | `acyl` | 1 | `()` | `parent_anchor_fields=("acyl_c_idx",)` |
| 3 | `azanide` | 4 | `()` | 氮负离子母体（P-72.2.2.2），类 4 阴离子，高于酸 |
| 4 | `cation` | 6 | `()` | 单核母体阳离子（P-73.1.1），高于酸 |
| 5 | `acid` | 7 | `(1,)` | 羧酸 |
| 6 | `oxoacid` | 9 | `(0,)` | 含氧酸中心 P/S/B 合一类，`oxo_kind` 由 L1 payload 归一 |
| 7 | `sulfonamide` | 11 | `(1,)` | P-41 类 11，与酰胺同组、排在酰胺之后 |
| 8 | `ester` | 9 | `()` | 环外附着取位次 |
| 9 | `acyl_halide` | 10 | `()` | |
| 10 | `amide` | 11 | `()` | 环外附着取位次 |
| 11 | `nitrile` | 14 | `()` | 环外附着取位次 |
| 12 | `aldehyde` | 15 | `()` | 环外附着取位次 |
| 13 | `ketone` | 16 | `()` | |
| 14 | `thione` | 16 | `(1,)` | C=S（P-64.6），与酮同类号，路径靠后让位于 C=O |
| 15 | `alcohol` | 17 | `(1,)` | 锚点 `surr_idx` |
| 16 | `thiol` | 17 | `(2,)` | 与醇同类号，路径靠后 |
| 17 | `amine` | 19 | `()` | 锚点 `surr_idx` |
| 18 | `heterane` | 36 | `()` | 非碳母体氢化物（P-41 类 21–39），低于胺 19、高于碳 40 |

类号越小越优先。`layer2/principal.PRINCIPAL_REGISTRY` 由 `FG_SPECS` 中 `p41` 非零的条目派生，经 `_spec_from_fg` 转成 `PrincipalFeatureSpec(PrincipalPriority(p41, path), expr, parent_anchor_fields)`。`PrincipalPriority` 是 `order=True` 的冻结 dataclass，字段为 `p41_class` 与 `p43_path`，因此可直接比较：先比类号，再按路径逐项比较。`select_principal_group` 用 `min(..., key=...)` 取出最小者，即优先级最高的类，再取该类全部 occurrence。

路径只在同类号内部起作用。`alcohol` 与 `thiol` 同为类 17，靠 `path` 分出先后：醇是 `(1,)`、硫醇是 `(2,)`，因此羟基优先作后缀，巯基让位。`ketone` 与 `thione` 同为类 16，酮是 `()`、硫代酮是 `(1,)`，C=O 优先作后缀、C=S 让位。`acid` 的 `p41=7` 与 `path=(1,)` 使它排在类 7 内的首位；`oxoacid` 的 `path=(0,)` 只在它留在类 9 时参与比较，升到类 7 后沿用的仍是同一路径（见下节）。`sulfonamide` 的 `path=(1,)` 让它排在 `amide`（`path=()`）之后，两者同属类 11。`heterane` 只占一个等级 36：P-41 类 21–39 整段低于胺、高于碳，元素间次序另由 L2 的 `constants.PARENT_SENIOR_ATOMS`（P-44.1.2）裁决。

`FG_SPECS` 的注册顺序即 L1 检测列表的产出顺序：`functional_group_inventory._FG_KEYS` 直接取 `tuple(sp.fg for sp in FG_SPECS)`，清单条目按此顺序展开。`_ANCHOR_KEYS` 同样由 `FG_SPECS` 过滤 `anchors` 非空的条目派生，只有声明了锚点的类别才在 occurrence payload 中收集锚点。`FunctionalGroupClass.NONE` 的值为 `alkane`，用于无主官能团的场合：`select_principal_group` 在 `eligible` 为空时返回 `PrincipalGroupSelection(FG.NONE, ())`，把母体交给烷烃处理。

`FgSpec.locant_source` 决定位次原子从哪里取：`attachment` 取 `principal_expression_facts` 中骨架内的附着原子，是默认值；`attachment_exocyclic` 只在环外表达时取，`ester`、`amide`、`nitrile`、`aldehyde` 用此值；`anchor_field` 取 `parent_anchor_fields` 首位语义字段，即固定 locant 1 的锚点。

## 类 7 的含氧酸漂移

`oxoacid` 在 `FG_SPECS` 中声明为 `p41=9`，但含氧酸的实际类别可随酸式漂移：同一个含氧酸中心，写成酸式或酸根（带 O⁻）时按 P-41 类 7 的「酸」参与竞争，写成中性酯时留在类 9。漂移后的等级常量是 `OXO_ACID_P41 = 8`，即类 7 内排在羧酸（`acid` 的 `7`）之后。

判定某个 occurrence 是否按酸式参与竞争由 `fg_registry.oxoacid_is_acid(payload)` 给出，三条判据依次为：

- 中心带 O⁻（`n_om > 0`）→ 按酸根处理，返回真。
- `oxo_kind` 落在 `OXO_ACID_KINDS`（当前为 `frozenset({"sulfonic"})`）→ 碳锚定且自带酸式氢，返回真。
- `oxo_kind` 落在 `OXO_ACID_KIND_BY_H`（当前为 `{"phosphonate": 1, "boronic": 1}`）→ 碳锚定 P/B 酸按氢数阈值判定，即 `n_oh >= 1` 时返回真。膦酸 / 硼酸与各自的酸氢酯因此落类 7，全酯化的膦酸酯 / 硼酸酯留在类 9。

漂移的消费点在 `layer2/principal._effective_priority(group_class, spec, inventory, mol=None)`。它是唯一允许类别漂移的地方，且只对 `FG.OXOACID` 生效。函数要求该类全部 occurrence 都满足 `oxoacid_is_acid`（`occurrences and all(...)`），即整类一致才升到 `PrincipalPriority(OXO_ACID_P41, spec.priority.p43_path)`；只要有一条不是酸式，整类保持 `spec.priority` 的类 9。路径 `p43_path` 原样带走，因此漂移只换类号、不改路径。`select_principal_group` 把 `_effective_priority` 作为 `min` 的 key，得到的就是漂移后的有效优先级。

四个 kind 的分类去向：

| `oxo_kind` | 判据 | 分类去向 |
| --- | --- | --- |
| `sulfonic` | 恒为酸式 | 类 7，等级 8，作主基团出磺酸后缀 |
| `phosphonate` | `n_oh >= 1` | 类 7，等级 8；`n_oh == 0` 时回类 9 |
| `boronic` | `n_oh >= 1` | 类 7，等级 8；`n_oh == 0` 时回类 9 |
| `phosphate` / `sulfate` | 带 O⁻ 时按酸根 | 中性时类 9；带 O⁻ 时类 7 |

`layer1/analyzer._OXO_CLASS_BY_KIND` 负责把 `oxo_kind` 映射到 P-41 类别键，`{"sulfonamide": "sulfonamide"}` 使磺酰胺从含氧酸合一类中分出去，落类 11，其余 kind 留在 `oxoacid` 类。

## 阴离子与阳离子的定级

`_effective_priority` 另有两处按电荷定级的漂移，与含氧酸漂移并列，都只在整类一致时生效。

- O⁻/S⁻ 型醇盐 / 硫醇盐：`FG.ALCOHOL` / `FG.THIOL` 的全部 occurrence 带 `anion_os` 标记（L1 `analyzer._mark_anion_os` 给中心为 O⁻/S⁻ 的条目置位），且 `principal._has_charged_acid` 判定分子内无阴离子型酸根时，整类升到 `PrincipalPriority(_ANION_OS_P41, ...)`，`_ANION_OS_P41 = 4`，即 P-41 表 4.1 类 4 阴离子高于酸（类 7）。若有羧酸根 / 含氧酸根在场，氧 / 硫负离子不抢母体，仍由酸作母体。
- 阳离子让位：`FG.CATION` 在 `inventory.has_anion` 为真时置底（`_CATION_ANION_GATED = 99`），不再作母体（类 4 > 类 6）。

`_has_charged_acid` 分别扫 `FG.ACID` 的 `surr_idx` 中是否有负电荷氧、`FG.OXOACID` 的 `n_om > 0`，命中即认为阴离子型酸在场。

L1 的 `functional_group_inventory.FunctionalGroupInventory.has_anion` 字段承载「分子内有负形式电荷原子」这一事实，`build_inventory` 由 `analyzer._has_negative_atom` 一次扫描后传入，避免重复遍扫。

## 三层接力

优先级数据从声明到成词穿过三层，每层只消费上一层给出的事实。

| 层 | 模块与函数 | 对优先级的消费方式 |
| --- | --- | --- |
| L1 声明 | `fg_registry.FG_SPECS` / `FgSpec.p41` / `p43_path` | 唯一声明处；同时提供 `oxoacid_is_acid` 判定酸式 |
| L1 仲裁 | `analyzer._arbitrate_parts` | 按 `p41` 比较，压制组合 FG、标记降级叶 |
| L2 选择 | `principal._effective_priority` / `select_principal_group` | 取有效优先级最小者作主基团 |
| L2 骨架 | `parent_skeleton.select_principal_skeletons` / `parent_select._reorder_p45_2` | P-44 筛选与 P-45.2 排序 |
| L5 出词尾 | `chain_engine._KIND_TABLE` / `_OXO_TAIL` | 按 kind 取后缀，含氧酸按 `oxo_kind` 与氢数查表 |

L5 的 `_KIND_TABLE` 以 kind 字符串为键，每个 `_Chain` 记录 `en_suf` / `zh_suf` 等词尾信息，其中 `fg` 字段指回 FG 键，把产物重新绑到 L1 的类别上：

| `_KIND_TABLE` 键 | `fg` 字段 | 词尾 |
| --- | --- | --- |
| `alcohol` | `alcohol` | `ol` / 醇 |
| `ketone` | `ketone` | `one` / 酮 |
| `thione` | `thione` | `thione` / 硫酮，`coda="ane"` |
| `alkane` | 无 | `ane` / 烷 |
| `acid` | 无 | `oic acid` / 酸 |
| `sulfonic` | `oxoacid` | `sulfonic acid` / 磺酸，`coda="ane"` |
| `sulfonate` | `oxoacid` | `sulfonate` / 磺酸 |
| `sulfonamide` | `sulfonamide` | 磺酰胺 |
| `sulfonyl_chloride` | `oxoacid` | 磺酰氯 |
| `ester` | `ester` | `oate` / 酸 |
| `phosphate` | `oxoacid` | `phosphate` / 磷酸 |
| `phosphonate` | `oxoacid` | `phosphonic acid` / 膦酸 |
| `sulfate` | `oxoacid` | `sulfate` / 硫酸 |
| `boronic` | `oxoacid` | `boronic acid` / 硼酸，`plain_hook=_oxoacid_tail` |
| `acyl` | `acyl` | `oyl` / 酰基 |
| `thiol` | `thiol` | `thiol` / 硫醇，`coda="ane"` |
| `amine` | `amine` | `amine` / 胺 |
| `aldehyde` | `aldehyde` | `al` / 醛 |
| `nitrile` | `nitrile` | `enitrile` / 腈 |
| `amide` | `amide` | `amide` / 酰胺 |
| `radical` | `radical` | `yl` / 基，`coda="an"` |

含氧酸中心母体的词尾不由 `en_suf` 直接给出，而由 `_OXO_TAIL` 查表：键是 `(oxo_kind, 中心酸式氢数)`，例如 `("phosphate", 3)` 取 `("phosphoric acid", "磷酸")`、`("phosphonate", 2)` 取 `("phosphonic acid", "膦酸")`、`("boronic", 2)` 取 `("boronic acid", "硼酸")`、`("sulfate", 1)` 取 `("hydrogen sulfate", "硫酸氢")`。表内覆盖 `phosphate` / `phosphonate` / `sulfate` / `boronic` 四族、氢数从满酸式递减到全取代共十三项：氢数等于该酸的最大酸式氢数时出「酸」义（`phosphoric acid` / 磷酸、`sulfuric acid` / 硫酸、`boronic acid` / 硼酸），氢数减少时逐级切到「二氢」「氢」义（`dihydrogen phosphate` / 磷酸二氢、`hydrogen phosphate` / 磷酸氢），氢数为 `0` 时出「酯」义（`phosphate` / 磷酸、`sulfate` / 硫酸、`boronate` / 硼酸）。查表函数 `_oxoacid_tail` 读 `numbered["parent"]` 的 `oxo_kind` 与 `n_oh`，表外返回 `None`，由调用方回落到常规 `en_suf`。

L1 的 `FunctionalGroupClass` 枚举与 `FG_SPECS` 的 `fg` 值一一对应，`sulfonamide` 是独立枚举值而非 `oxoacid` 的子类，因此 L5 侧磺酰胺走 `fg="sulfonamide"` 的分支，与酰胺同组而不与含氧酸同组。`_KIND_TABLE` 中 `fg` 字段缺失的 kind（`alkane`、`acid`）表示该 kind 不需要回绑 FG 键，`acid` 的词尾由 `variant` 中的保留名（如 `oxalic acid` / 草酸）承担。

## 压制与降级

L1 内部仲裁在 `analyzer._arbitrate_parts(parts, mol, has_anion)` 完成，涉及三个常量。进入压制判定前先做两道同族同位次去重：`_drop_mixed_anion_os` 在同一 `alcohol` / `thiol` 族内同时存在中性羟基与 O⁻/S⁻ 时只留负离子条目（酚盐 / 硫醇盐优于中性酚 / 硫醇，P-66.1.1.4）；`_drop_mixed_anion_acids` 在同一 `acid` / `oxoacid` 族内只留带负电荷氧的条目（阴离子优先于中性酸，P-65.6.2.3.1）。`has_anion` 由 `_has_negative_atom` 一次扫描后传入，为真时把 `cation` 移出存在性集合（类 4 阴离子 > 类 6 阳离子，避免阳离子误压制酸）。

`_SUPPRESSIBLE` 是可被更高优先级 FG 整体压制的组合 FG 集合：`{"acid", "ester", "acyl_halide", "amide", "nitrile", "aldehyde"}`，即组合羰基加腈。仲裁先由 `FG_SPECS` 造出 `p41` 字典，再收集「存在」的 FG 集合；对 `_SUPPRESSIBLE` 中每个类别，只要存在另一个类号更小的 FG，该类就退出。退出方式分两种：`_LEAF_DEMOTED = ("acid", "nitrile")` 中的酸与腈整组碳被排除出主链，其 occurrence 标 `demoted=True`，成为 P-61.1.3 的前缀叶（`carboxy` / `cyano`）；其余组合 FG（酯/酰胺/醛/酰卤）的条目直接清空为「氧代」，羰基碳留在链内，只把 O 作 `oxo` / `formyl` 前缀，由 L3 锚定叶识别。

`_PRESENCE_SKIP = {"oxoacid", "sulfonamide"}` 把含氧酸与磺酰胺排除在存在性判定之外。这两个类若计入 `present`，会参与对 `_SUPPRESSIBLE` 的比较并改写压制结果，因此仲裁阶段不看它们。

降级状态随清单传递：`functional_group_inventory.FunctionalGroupOccurrence.demoted` 记录该标记，`FunctionalGroupInventory.occurrences` 与 `demoted_entries` 分别返回未降级与已降级的条目。`select_principal_group` 只用 `not entry.demoted` 的条目作候选，`_chain_candidates` 则通过 `_demoted_leaf_carbons` 把降级叶的碳排除出开链枚举。

`build_inventory(lists, mol, demoted, has_anion)` 接收 `_arbitrate_parts` 给出的 `demoted` id 集（形如 `f"{fg}:{i}"`），组装 occurrence 时按 id 命中即置位。因此降级状态在 L1 出口就已固化进清单，L2 与后续各层只读取该结果，不重算。

三类状态因此可以并存于同一分子：主基团（`select_principal_group` 选出的类，出后缀）、前缀 FG（未压制但也未当选的类，出前缀）、降级叶（被压制且标 `demoted` 的酸/腈，出 `carboxy` / `cyano` 前缀）。`_SUPPRESSIBLE` 中酯/酰胺/醛/酰卤被压制时条目清空、不标 `demoted`，其氧化态由 L3 锚定叶从骨架中重新识别，所以它们的碳留在主链内参与编号。

## P-45.2 排序与缩合磷酸次序

P-44 拓扑规则只活在 `parent_skeleton.select_principal_skeletons` 的筛选谓词里。该函数按顺序调用：`keep_max_principal_coverage` 取主基团覆盖度最大者，`keep_p44_1_2` 在混合拓扑时取 senior 元素，随后按拓扑分流——纯开链走 `keep_p44_3`（杂原子数、原子数、元素计数取最大），含环走 `keep_p44_2`，最后统一由 `keep_p44_4_unsaturation` 按不饱和度键收尾。

P-45.2 排序在 `parent_select._reorder_p45_2`。它对候选列表做稳定排序，排序键是三元组，比较方向为：

1. `_p45_2_prefix_count(info, cand)` 降序——该候选 `owned_atoms` 边界外的 claim 个数，即前缀取代基团数目（P-45.2.1），少者优先。
2. `_condensed_rank(info, cand)` 降序——缩合含氧酸（磷酸 / 硫酸）的链内桥氧数（P-67.2.1），多者优先。
3. 原始序升序——保持枚举顺序，作为平局兜底。

`_condensed_rank` 只对 `oxo_kind ∈ ("phosphate", "sulfate")` 的候选计算，其余候选恒返回 `0`，因此在第一键平局时缩合磷酸 / 硫酸候选才能凭第二键胜出。它从候选的 `covered_principal_ids` 与 `principal_occurrences` 中取出 `oxo_z`，对每个中心调用 `analyzer._oxo_bridge_arms` 并取最大值。`_oxo_bridge_arms` 统计中心原子的 O-桥氧数：中心须为 P 或 S（原子序数 15 / 16，否则恒返回 0），逐个数非双键氧邻居中是否还连着同种中心原子，是则计一个桥氧。链内桥氧越多，说明该中心越是多核（P-O-P / S-O-S）链的功能母体，排序越靠前。`_reorder_p45_2` 在 `tied=True` 时只保留第一键取最大值的候选，供 `select_parent` 返回并列最优。

缩合含氧酸的候选资格在 L1 阶段就已收窄：`analyzer.oxoacid_entries` 在识别含氧酸中心后做两道过滤。第一道要求链上留有全酸式末端——若没有任何中心的 `n_oh + n_om >= 2`，带 P-O-P / S-O-S 桥氧的中心被剔除，即全酯化的缩合链不按功能母体识别（P-67.2.1）；此时整链也不退为取代基，至少留酸式氧（`n_oh + n_om`）最多的链节作母体。第二道取全部中心里最大的 `n_om`，把 `n_om` 更小且带桥氧的链内中心剔除，让位给质子化更少的酸中心，避免酸根被写成前缀。只有通过这两道的候选，才会带着 `oxo_z` 进入 `_condensed_rank` 的桥氧计数。

`select_parent` 的调用链是 `_collect_candidates` → `_finalize_ranked` → `_reorder_p45_2(tied=True)` → `_reorder_oxo_ester_side`：先由规则驱动收集候选并去重，再由 `_finalize_ranked` 补齐契约、词干与编号并固化 `owned_atoms`，接着交 P-45.2 排序取并列最优，最后对缩合磷酸酯的并列候选做母体侧定向——`_oxo_ester_side_score` 按（臂上芳香含氮环数, -臂内氧数）降序取「糖 / 多元醇 / 甘油」侧作母体（P-67.1.3），仅 `oxo_kind == "phosphate"` 的候选参与。`_p45_2_prefix_count` 依赖的 `owned_atoms` 正是在 `finalize_parent_ownership` 中固化，因此排序键所用的 claim 计数与最终归属一致。

## 相关文件

- `src/namepredict/layer1/fg_registry.py`：`FgSpec`、`FG_SPECS`、`oxoacid_is_acid`、`OXO_ACID_P41`、`OXO_ACID_KINDS`、`OXO_ACID_KIND_BY_H`
- `src/namepredict/layer1/functional_group_inventory.py`：`FunctionalGroupClass`、`FunctionalGroupOccurrence`、`FunctionalGroupInventory`、`OXO_FG_CLASSES`、`FG_ATOM_FNS`、`_FG_KEYS`、`_ANCHOR_KEYS`
- `src/namepredict/layer1/analyzer.py`：`_arbitrate_parts`、`_SUPPRESSIBLE`、`_PRESENCE_SKIP`、`_LEAF_DEMOTED`、`_OXO_CLASS_BY_KIND`、`_oxo_bridge_arms`、`_has_negative_atom`、`_drop_mixed_anion_os`、`_drop_mixed_anion_acids`、`_mark_anion_os`、`_OXO_KIND_P`、`_OXO_KIND_S_ARM`、`boronic_entries`
- `src/namepredict/layer2/principal.py`：`PrincipalPriority`、`PRINCIPAL_REGISTRY`、`_effective_priority`、`_has_charged_acid`、`_ANION_OS_P41`、`_CATION_ANION_GATED`、`select_principal_group`
- `src/namepredict/layer2/parent_skeleton.py`：`select_principal_skeletons`、`keep_p44_1_2`、`keep_p44_2`、`keep_p44_3`、`keep_p44_4_unsaturation`
- `src/namepredict/layer2/parent_select.py`：`_p45_2_prefix_count`、`_condensed_rank`、`_reorder_p45_2`、`_reorder_oxo_ester_side`、`_oxo_ester_side_score`、`select_parent`
- `src/namepredict/layer5/chain_engine.py`：`_KIND_TABLE`、`_OXO_TAIL`、`_oxoacid_tail`

相关页面：[[architecture/layer1-analyzer]]、[[architecture/layer2-parent-selector]]、[[concepts/atom-ownership]]、[[reference/core-data-contracts]]、[[guides/adding-new-functional-group]]。
