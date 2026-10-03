# 指南：新增官能团 (Adding a Functional Group)

> 扩展落地清单 | 关联：[[architecture/overview]]、[[concepts/functional-group-priority]]、[[concepts/atom-ownership]]、[[reference/core-data-contracts]]

主线为「检测 → 归一 → 注册 → 母体化 → 命名」五段。主示例取碳锚定的磺酸族 `sulfonic`：中心 S 不是碳骨架成员，经一条直连碳臂接入母体。文末给中心自任母体的对照示例。

| 段落 | 落点 | 文件 |
|---|---|---|
| 检测 | `FG_SMARTS` | `layer1/fg_local_smarts.py` |
| 归一 | `_oxo_kind` / `_oxoacid_entry` | `layer1/analyzer.py` |
| 分类 | `_OXO_CLASS_BY_KIND` | `layer1/analyzer.py` |
| 注册 | `FG_SPECS` / `FunctionalGroupClass` / `FG_ATOM_FNS` | `layer1/fg_registry.py`、`layer1/functional_group_inventory.py` |
| 母体化 | `_kind_fg_atoms` / `OXO_FG_CLASSES` / `OXO_CENTER_KINDS` | `layer2/parent_select.py`、`layer1/functional_group_inventory.py`、`constants.py` |
| 命名 | `_KIND_TABLE` / `_OXO_TAIL` / `_BENZENE_RETAINED` | `layer5/chain_engine.py` |

## 1 局部检测：`FG_SMARTS`

`FG_SMARTS` 每条为 `(FG 键, SMARTS)`，**模式首原子即该 FG 的中心原子**；`match_local_fg` 返回 `{FG 键: [原子元组, ...]}`，按中心原子升序、同名多条取并集。磺酸族的表项（S 带两个 =O，另两臂为酸式氧与直连碳臂）：

```python
("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS + _C_ARM)
```

臂词常量按需复用（不另起一份）：`_ACID_O`、`_NOT_ACID`、`_NOT_ACYCLIC_ESTER`、`_NOT_HALO`、`_ONE_C`、`_RING_HET`、`_O_PHOS`、`_O_ARM`、`_O_S_ARM`、`_HALO_ARM`、`_N_ARM`、`_C_ARM`、`_SNY_ARM`。同族按臂型各占一条：`_HALO_ARM` → 磺酰卤、`_N_ARM` → 磺酰胺、`_O_ARM` → 磺酸酯、`_O_PHOS + _O_ARM` → 硫酸酯、`_O_PHOS * 2` → 硫酸（含硫酸根）、`_O_PHOS + _O_S_ARM` / `_O_S_ARM * 2` → 多硫酸链端节 / 中节（`_O_S_ARM` 为 S-O-S 桥臂）。同一个 FG 键也可由多条 SMARTS 覆盖不同双键杂原子：`amide` 键有 C=O / C=S / C=N 三条（`=O` 酰胺 / 硫代酰胺 / 脒同属 P-41 类 11），歧义由 L2 写入的 `parent["amide_z"]` 消解、L5 按它切换词尾。`_compile` 按键归并编译，无法解析即报错。

## 2 非局部归一：`analyzer` 的含氧酸管线

需要看整分子环境的 FG 才走本段。含氧酸中心（B/P/S 非碳骨架成键）共用一条管线：

1. `_oxo_partition` 数中心邻居角色（双键氧 / 羟基氧 / 阴离子氧 / O-臂 / 碳 / 卤素 / 氮 / 硫），返回 `oxo` / `oh` / `om` / `o_arm` / `o_arm_root` / `c` / `hal` / `n` / `s` 分桶；
2. `_oxo_kind` 归一 `oxo_kind`：P 走 `_OXO_KIND_P`（键为 `(双键氧数, 有无直连碳臂)`），S 走 `_OXO_KIND_S_ARM`（键为另一臂角色 `halo` / `n` / `acid` / `o_arm`，值即 kind）；S 无直连碳臂时不查该表，`n_oh + n_om + o_arm == 2` 即 `sulfate`（硫酸、硫酸根、多硫酸链节）；B 不登记 `FG_SMARTS`，由 `boronic_entries` 按元素扫描，恰一个碳臂 + 两个酸式氧即 `boronic`；
3. `_oxoacid_entry` 组装条目，入口 `oxoacid_entries`（B 另由 `boronic_entries` 补入）。

两条分支决定条目形态：

| 分支 | 判据 | 条目键 |
|---|---|---|
| 中心自任母体 | `oxo_kind in OXO_CENTER_KINDS` | `oxo_z` / `center_idx`（= 中心）/ `oxo_kind` / `n_oh` / `n_om` |
| 碳锚定 | 直连碳臂唯一（`len(c_roots) == 1`） | 同上，但 `center_idx` 取该碳，另带 `surr_idx`（只收阴离子氧，供 L2 补 anion 标志） |

中心自任母体另需 `_arm_single_attach` 单点回接与整分子纯度校验（重原子 = core ∪ 臂）。

`_OXO_KIND_S_ARM`：S 带直连碳臂时，另一臂角色 → kind。

| 臂角色键 | 判据 | kind |
|---|---|---|
| `halo` | 卤素臂 | `sulfonyl_chloride` |
| `n` | 非环氮臂 | `sulfonamide` |
| `acid` | 酸式氧（OH / O⁻） | `sulfonic` |
| `o_arm` | O-R 臂 | `sulfonate` |

加同族成员时只在本表补一行，不新写检测器。`_oxo_bridge_arms` 数中心 P/S 的 O-桥氧。`oxoacid_entries` 另做两项收尾：缩合含氧酸（P-O-P / S-O-S）须链上留有全酸式末端才按功能母体识别；无全酸式末端时至少留酸式氧最多的链节作母体（否则母体落到甲烷）；链内中心让位于更少质子化的中心。

## 3 P-41 分类：`_OXO_CLASS_BY_KIND`

`oxoacid_lists` 用 `_OXO_CLASS_BY_KIND` 把条目归位到 FG 键：键存在即落该键，否则落 `oxoacid`。磺酰胺经 `_OXO_KIND_S_ARM["n"]` 得 kind `sulfonamide`，再由本表落 `"sulfonamide"`（与酰胺同属类 11）；新 kind 落别的类时改这张表。

## 4 注册元数据

`layer1/fg_registry.FG_SPECS` 加一条 `FgSpec`（字段按源码共 7 个）：

| 字段 | 填什么 |
|---|---|
| `fg` | FG 键 / 类别值，同时是 L1 条目键与 L4 位次 kind |
| `p41` | P-41 主官能团等级，`0` = 非主官能团；越小越优先 |
| `path` | P-43 同类内路径，与 `p41` 合成 `PrincipalPriority` |
| `expr` | 表达类型，取值见 `PrincipalExpression` |
| `anchors` | payload 中承载母体连接原子的键；空 = 不收集 |
| `parent_anchor_fields` | `(单, 复)` 锚点字段名，仅固定 locant 1 的扁平字段 |
| `locant_source` | `attachment` / `attachment_exocyclic` / `anchor_field` |

磺酸族与磺酰胺各一条：

```python
FgSpec("oxoacid", p41=9, path=(0,), anchors=("center_idx",))
FgSpec("sulfonamide", p41=11, path=(1,), anchors=("center_idx",))
```

`layer1/functional_group_inventory.FunctionalGroupClass` 加成员，值必须与 `FgSpec.fg`、`FG_SMARTS` 键逐字一致——`_ANCHOR_KEYS` 与 `_one` 都做 `FunctionalGroupClass(key)`，缺成员即失败。

`FG_ATOM_FNS` 只在特征原子不走「`center_idx` ∪ `surr_idx`」通用形态时登记：`oxoacid` 与 `sulfonamide` 共用 `_oxoacid_atoms`，即「锚点碳 + 中心 B/P/S + 中心的非碳非氢邻居」；`azanide`、`heterane` 与 `cation` 共用 `_cation_atoms`（只占中心本身），`nitrile` 用 `_nitrile_atoms`（腈碳 + 三键氮，R 侧不留）。

## 5 P-41 酸式漂移（可选）

该类可整体升为 P-41 类 7 竞争时，用 `fg_registry.oxoacid_is_acid` 判定：中心带 O⁻、或 `oxo_kind in OXO_ACID_KINDS`、或 `n_oh >= OXO_ACID_KIND_BY_H[oxo_kind]`（当前含 `phosphonate` / `boronic`）；命中时 `layer2/principal._effective_priority` 换成 `PrincipalPriority(OXO_ACID_P41, path)`（`OXO_ACID_P41 = 8`，排在羧酸之后）。三张表同在 `fg_registry`，按新 kind 补行即可。电荷型定级同在此函数：醇/硫醇盐经 `anion_os` 升类 4（`_ANION_OS_P41`），阳离子在 `inventory.has_anion` 时置底（`_CATION_ANION_GATED`）。

## 6 母体化接线（L2）

注册完成后 L2 自动接线，只有例外情形才动手：

- `select_principal_group` 取 `PrincipalPriority` 最小的类别作主基团，低优先级 FG 一律成为取代基（互斥由「只选单个最高优先级类」实现）。
- `express_chain_principal` / `express_ring_principal` 按骨架拓扑表达主基团；`_chain_kind` 给出 kind，含氧酸类取 `_oxo_kind_of`，要求该类全部 occurrence 的 `oxo_kind` 一致，否则该类不支持。
- `_chain_oxoacid_fields` 写入 `oxo_kind` / `n_oh` / `n_om`；`oxo_kind in OXO_CENTER_KINDS` 时另做盐门控（金属数与 `n_om` 匹配，不过则返回 `None` 使候选作废），碳锚定类无门控。
- `_facts` 汇总 `group_class` / `multiplicity` / `relation`（`in_skeleton` / `exocyclic`）/ `occurrence_ids` / `characteristic_atoms` / `anchor_atoms` / `attachment_atoms`；多重度由 occurrence 个数承载，不要另造 kind。
- `assemble` 依次接 hydro 前缀、环内阳离子后缀、核素描述符（`join_isotope_descriptor`）、`join_kind_name`、环内 λ 杂原子 =O 的 oxide 分离、硫族/腈别名切换、阴离子、E/Z、R/S；盐后缀由 `namer._apply_salt_suffix` 收尾。

例外挂钩点：

| 需求 | 挂钩点 |
|---|---|
| 位次由固定 locant 1 的扁平字段承载 | `_SEMANTIC_ANCHOR_FGS` / `_semantic_anchor_fields` 放行并写字段，配 `FgSpec.parent_anchor_fields` |
| L5 整名需要新的计数或盐元数据 | `express_chain_principal` 中按 `kind` 加 producer |
| 该 FG 被更高优先级类压制 | `analyzer._SUPPRESSIBLE`（组合 FG）/ `_LEAF_DEMOTED`（降为前缀叶） |

## 7 命名

`layer5/chain_engine._KIND_TABLE` 注册该 kind 的 `_Chain` 条目。字段按源码列全：

```text
词形 kind / en_suf / zh_suf / coda
位次 fg / need / no_loc / omit_rule / yl_loc_omit / zh_loc_omit / ene_loc_omit / yne_loc_omit
不饱和 ene_seg / yne_seg / ene_base / yne_suf / unsat_polyol
形态 ez_ene / ez_ene_multi / wrap / cyclic / cyclic_unsat / aromatic / stem / zh_full
保留名 plain_maps / plain_fn / plain_hook / variant
```

磺酸的条目：

```python
"sulfonic": _Chain(kind="sulfonic", en_suf="sulfonic acid", zh_suf="磺酸", coda="ane",
                   fg="oxoacid", need=1, omit_rule=_omit_term_locant,
                   ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent)
```

- `fg` 是**位次记录的 kind**（= `FgSpec.fg`），不一定等于 `kind`：`sulfonic` 的 `fg` 填 `"oxoacid"`（L1 条目键与 L4 记录键都是 `oxoacid`），`sulfonamide` 两者同名。
- `need=1` + `omit_rule=_omit_term_locant` 让单取代短链省位次；多基后缀由 `_generated_mult_fields` 按 `facts.multiplicity` 自动生成，不要为多重度另造 kind。

载荷字段的取值含义：

| 字段 | 取值 / 作用 |
|---|---|
| `coda` | 饱和词干与后缀间的词尾，`"ane"` 配元音开头后缀时省 `e` |
| `no_loc` | `"plain"` 出普通名；`"none"` 位次记录缺失即返回 `None` |
| `omit_rule` | `(n, loc, omit) -> bool`，`True` 即省略位次 |
| `plain_hook` | `(numbered) -> pair`，词尾由分子计数决定，不查碳数词表 |
| `variant` | `{scaffold_id: {multiplicity: 覆盖字段}}`，键 `None` 覆盖无 scaffold 情形 |
| `plain_maps` / `plain_fn` | 俗名表 / 派生命名；皆空时回落词干拼接 |
| `ene_base` / `yne_suf` | 融合式烯/炔基座；`ene_seg` / `yne_seg` 为段式段形 |
| `unsat_polyol` | 多 FG 词干支持烯/炔插入（多醇 / 多硫醇） |
| `cyclic` / `cyclic_unsat` | 恒加环前缀 / 仅在不饱和段加环 |

其余词尾落点：

| 场景 | 落点 |
|---|---|
| 中心自任母体的词尾 | `_OXO_TAIL`，键为 `(oxo_kind, 中心酸式氢数)`，由 `plain_hook=_oxoacid_tail` 消费；碳锚定族不查此表 |
| 苯环单取代保留名 | `_BENZENE_RETAINED`（`sulfonic` → benzenesulfonic acid / 苯磺酸） |
| 其他 scaffold 单取代保留名 | `_Chain.variant` 的 `{scaffold_id: {multiplicity: 覆盖字段}}` |
| 酰胺族按双键杂原子换词尾 | `parent["amide_z"]` 记 `=O`/`=S`/`=N` 的杂原子序数（`principal_expression` 写入）；L5 据它把 `amide` 词尾换成 `thioamide` / 硫代酰胺 或 `imidamide` / 亚氨酰胺，环外经 `_exo_ring_spec` 出 `carbothioamide` / `carboximidamide` |
| 环外主基系统名 | `constants.EXO_RING_SUF` 加 `group_class → (单, 复)`；由 `_exo_ring_spec` 改写词尾/词干/位次 |
| 单碳母体带两个杂原子 / 官能团碳直连 N 的词尾 | `assembler._c1_retained`（`urea` / `thiourea` / `guanidine` / `carbonate`）与 `_c1_amino`（`carbamate` / `carbamic acid` / `carbamoyl`）；`locant_kind` 记保留名，`subs_consumed` 表示取代基已并入母体名 |
| 作取代基前缀时的保留名 | `src/namepredict/tools/anchored_table._REGISTRY` 加 `RetainedSubstituent`（如 `sulfo` / `sulfonato` / `borono` / `trimethylsilyl` / `selanyl` / `thiocyanato` / `diazonio`），`anchored` 用带 `*` 的 canonical SMILES；`formyl` 的中文名为「甲酰基」。`-C(=O)-N(R)(R')` 不走注册表，由 `carbamoyl_prefix_name` 收成 `<N-取代基>carbamoyl` 或环胺的 `<母体>-<位次>-carbonyl`（P-65.2.1.5） |

- 保留名母体的位次口径：`constants.N_LOCANT_KINDS` 的 kind 让 `assembler_prefixes._omit_sub_locants` 判「位次不可省」（`1,3-二甲基脲`）；`locant_kind == "carbamic_acid"` 反之（`dimethylcarbamic acid`）。

## 8 原子归属

| 条件 | 落点 |
|---|---|
| 特征原子跨两跳，或中心不是碳骨架成员 | `functional_group_inventory.OXO_FG_CLASSES` 加该 FG 类别；`parent_select._kind_fg_atoms` 命中后按 `covered_principal_ids` 整组纳入特征原子 |
| 中心自任母体、臂应退为取代基 | `constants.OXO_CENTER_KINDS` 加该 kind（`phosphate` / `phosphonate` / `sulfate` / `boronic`）；`finalize_parent_ownership` 把 `owned_atoms` 取成纯 FG 原子，链仅供编号 |

碳锚定的 `sulfonic` 不入 `OXO_CENTER_KINDS`。

## 9 自动跟进的派生表

以下由注册表投影，无需手改：

- `layer1/functional_group_inventory` 的清单键 `_FG_KEYS` 与 `_ANCHOR_KEYS` ← `FG_SPECS`；
- `layer2/principal.PRINCIPAL_REGISTRY` ← `p41 != 0` 的 `FgSpec`；
- L4 `layer4/locant_calc._FG_LOCANTS` ← `FG_SPECS`；
- 母体 kind 由 `principal_expression._chain_kind` 按 FG 类别现算，含氧酸类取自 L1 的 `oxo_kind`。

## 改动点 → 文件 → 是否派生

| 改动点 | 文件 | 派生自动 |
|---|---|---|
| `FG_SMARTS` | `layer1/fg_local_smarts.py` | 手动 |
| `_OXO_KIND_S_ARM` / `_OXO_CLASS_BY_KIND` | `layer1/analyzer.py` | 手动（仅含氧酸族） |
| `FG_SPECS` | `layer1/fg_registry.py` | 手动 |
| `FunctionalGroupClass` / `FG_ATOM_FNS` | `layer1/functional_group_inventory.py` | 手动 |
| 清单键 / `_ANCHOR_KEYS` | `layer1/functional_group_inventory.py` | 派生 |
| `_kind_fg_atoms`（`OXO_FG_CLASSES`） | `layer2/parent_select.py`、`layer1/functional_group_inventory.py` | 手动（仅跨跳特征原子） |
| `PRINCIPAL_REGISTRY` / kind | `layer2/principal.py`、`principal_expression.py` | 派生 |
| `OXO_CENTER_KINDS` / `EXO_RING_SUF` | `constants.py` | 手动（仅中心母体 / 环外表达） |
| `_KIND_TABLE` / `_OXO_TAIL` / `_BENZENE_RETAINED` | `layer5/chain_engine.py` | 手动 |
| `_c1_retained` / `_c1_amino` | `layer5/assembler.py` | 手动（仅 C1 双杂原子 / 官能团碳直连 N） |
| `parent["amide_z"]` 词尾切换 | `layer2/principal_expression.py`（写入）、`layer5/assembler.py` / `chain_engine._exo_ring_spec`（消费） | 手动（仅酰胺族双键杂原子置换） |
| `carbamoyl_prefix_name` | `tools/anchored_table.py` | 手动（仅 `-C(=O)-N(R)(R')`） |
| `_FG_LOCANTS` | `layer4/locant_calc.py` | 派生 |
| `_FG_GROUP` | `layer4/locant_calc.py` | 手动（仅位次省略特例） |
| `_REGISTRY` | `src/namepredict/tools/anchored_table.py` | 手动（仅作取代基前缀时） |

## 常见陷阱

1. **漏注册 `FgSpec`**：清单键集与 `_ANCHOR_KEYS` 都来自 `FG_SPECS`，缺该类则 inventory 不可见，下游整链断掉。
2. **`p41` 填错**：主基团胜负被静默改写。同类内用 `path` 排先后（`sulfonamide` 的 `path=(1,)` 排在 `amide` 之后）。
3. **特征原子范围决定 L3 提取边界**：范围过宽会抢走本该作取代基的碳；`_oxoacid_atoms` 只收中心的非碳非氢邻居与锚点碳，碳臂留给链/取代基侧。
4. **`_KIND_TABLE` 与 `_OXO_TAIL` 分离**：前者给基底词形，后者给中心母体按酸式氢数查的词尾。
5. **`_Chain.fg` 不等于 `_Chain.kind`**：`fg` 必须等于 `FgSpec.fg`（位次记录 kind），`sulfonic` 即反例。
6. **环形 / 无环、取代式 / 保留名的分叉要同时覆盖**：保留名走 `variant` / `_BENZENE_RETAINED`，环外走 `EXO_RING_SUF`，位次省略走 `_FG_GROUP`。
7. **派生表不要手写**：清单键、`_ANCHOR_KEYS`、`PRINCIPAL_REGISTRY`、kind、`_FG_LOCANTS` 都由注册表投影，另立一份必脱节。

## 验证

1. 简单 FG 作母体：最短链、不同链长各一条，断言词干与位次省略；
2. 不饱和链：验证段式引擎的烯/炔段形态；
3. 与更高优先级 FG 共存时退为取代基前缀，与其他取代基叠加时验证编号与位次；
4. 环骨架与环外表达各一条，确认保留名与 `EXO_RING_SUF` 分支被走到；
5. C1 双杂原子保留名（脲 / 硫脲 / 胍 / 碳酸酯）与 `-C(=O)-N(R)(R')` 取代基各一条。

双语（en / zh）都要断言，用例落 `tests/unit/`；pytest 并行运行，`test_architecture_contracts.py` 校验 `_detect_parts` 键集不越出 `FG_SPECS` 的 `fg` 集。

## 对照示例：中心自任母体（`phosphate`）

与碳锚定路径的三点差异：

- `_oxoacid_entry` 走 `OXO_CENTER_KINDS` 分支：`center_idx` 即中心原子，臂须单点回接并通过整分子纯度校验，条目不带 `surr_idx`；
- `owned_atoms` 只取 FG 特征原子（`finalize_parent_ownership` 依 `OXO_CENTER_KINDS` 判定），链仅供编号，臂退为取代基；
- 词尾走 `plain_hook=_oxoacid_tail` + `_OXO_TAIL`，按 `(oxo_kind, n_oh)` 出 phosphoric acid / dihydrogen phosphate 等，O 侧臂由 `join_kind_name` 按 `ESTER_O_SIDE_KINDS` 拼接。

中心 kind 归入 `OXO_CENTER_KINDS` 后，`namer._apply_salt_suffix` 自动跳过通用盐后缀，盐形态由该类自行组装。

相关页面：[[architecture/layer1-analyzer]]、[[architecture/layer2-parent-selector]]、[[architecture/layer5-name-assembly]]、[[concepts/functional-group-priority]]、[[guides/adding-new-ring-system]]。
