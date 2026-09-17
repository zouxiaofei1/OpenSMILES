# 指南：新增官能团 (Adding a Functional Group)

> 扩展落地清单 | 关联：[[architecture/overview]]、[[concepts/functional-group-priority]]、[[concepts/atom-ownership]]、[[reference/core-data-contracts]]

主线为「检测 → 归一 → 注册 → 母体化 → 命名」五段。主示例取碳锚定的磺酸族 `sulfonic`：中心 S 不是碳骨架成员，但经一条直连碳臂接入母体，属碳锚定路径。文末给中心自任母体的对照示例。

| 段落 | 落点 | 文件 |
|---|---|---|
| 检测 | `FG_SMARTS` | `layer1/fg_local_smarts.py` |
| 归一 | `_oxo_kind` / `_oxoacid_entry` | `layer1/analyzer.py` |
| 分类 | `_OXO_CLASS_BY_KIND` | `layer1/analyzer.py` |
| 注册 | `FG_SPECS` / `FunctionalGroupClass` / `FG_ATOM_FNS` | `layer1/fg_registry.py`、`layer1/functional_group_inventory.py` |
| 母体化 | `_WHOLE_FG_ATOMS` / `OXO_CENTER_KINDS` | `layer2/parent_select.py`、`constants.py` |
| 命名 | `_KIND_TABLE` / `_OXO_TAIL` / `_BENZENE_RETAINED` | `layer5/chain_engine.py` |

## 1 局部检测：`FG_SMARTS`

`FG_SMARTS` 是 `tuple[tuple[str, str], ...]`，每条 `(FG 键, SMARTS)`，**模式首原子即该 FG 的中心原子**。`match_local_fg` 统一返回 `{FG 键: [原子元组, ...]}`，按中心原子升序、同名多条取并集。磺酸族的表项（S 带两个 =O，另两臂为酸式氧与直连碳臂）：

```python
("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS + _C_ARM)
```

臂词常量按需复用，不另起一份：`_ACID_O`、`_NOT_ACID`、`_NOT_ACYCLIC_ESTER`、`_NOT_HALO`、`_ONE_C`、`_RING_HET`、`_O_PHOS`、`_O_ARM`、`_HALO_ARM`、`_N_ARM`、`_C_ARM`、`_SNY_ARM`。同族按臂型各占一条：`_HALO_ARM` → 磺酰卤、`_N_ARM` → 磺酰胺、`_O_ARM` → 磺酸酯、`_O_PHOS + _O_ARM` → 硫酸酯。`_compile` 按键归并编译，无法解析的模式立即报错。

## 2 非局部归一：`analyzer` 的含氧酸管线

需要看整分子环境的 FG 才走本段；不看环境的 FG 跳过。含氧酸中心（P/S 非碳骨架成键）共用一条管线：

1. `_oxo_roles` 数中心邻居角色（双键氧 / 羟基氧 / 阴离子氧 / O-臂 / 碳 / 卤素 / 氮 / 硫）；
2. `_oxo_kind` 查表归一 `oxo_kind`：P 走 `_OXO_KIND_P`（键为 `(双键氧数, 有无直连碳臂)`），S 走 `_OXO_KIND_S_ARM`（键为另一臂角色 `halo` / `n` / `acid` / `o_arm`，值即 kind）；
3. `_oxoacid_entry` 组装条目，入口 `oxoacid_entries`。

两条分支决定条目形态：

| 分支 | 判据 | 条目键 |
|---|---|---|
| 中心自任母体 | `oxo_kind in _OXO_Z_ANCHORED` | `oxo_z` / `center_idx`（= 中心）/ `oxo_kind` / `n_oh` / `n_om` |
| 碳锚定 | 直连碳臂唯一（`len(c_roots) == 1`） | 同上，但 `center_idx` 取该碳，另带 `surr_idx`（只收阴离子氧，供 L2 补 anion 标志） |

中心自任母体另需 `_arm_single_attach` 单点回接校验与整分子纯度校验（重原子 = core ∪ 臂，排除臂间成环、焦磷酸）。

`_OXO_KIND_S_ARM` 是「S 带直连碳臂时，另一臂角色 → kind」的对照表：

| 臂角色键 | 判据 | kind |
|---|---|---|
| `halo` | 有卤素臂 | `sulfonyl_chloride` |
| `n` | 有非环氮臂 | `sulfonamide` |
| `acid` | 有酸式氧（OH / O⁻） | `sulfonic` |
| `o_arm` | 有 O-R 臂 | `sulfonate` |

加同族成员时只在本表补一行，不新写检测器。`oxoacid_entries` 另做两项收尾：缩合磷酸（P-O-P）须链上留有全酸式末端才按功能母体识别；链内 P 让位于更少质子化的中心。

## 3 P-41 分类：`_OXO_CLASS_BY_KIND`

`oxoacid_lists` 用 `_OXO_CLASS_BY_KIND` 把条目归位到 FG 键：键存在即落该键，否则落 `oxoacid`。磺酰胺经 `_OXO_KIND_S_ARM["n"]` 得 kind `sulfonamide`，再由本表落 `"sulfonamide"`（与酰胺同属类 11）。新 kind 需要落别的类时改这张表。

## 4 注册元数据

`layer1/fg_registry.FG_SPECS` 加一条 `FgSpec`，字段按源码共 7 个：

| 字段 | 填什么 |
|---|---|
| `fg` | FG 键 / 类别值，同时是 L1 条目键与 L4 位次记录的 kind |
| `p41` | P-41 主官能团等级，`0` = 非主官能团；越小越优先 |
| `path` | P-43 同类内路径，与 `p41` 合成 `PrincipalPriority` |
| `expr` | 表达类型，取值见 `PrincipalExpression` |
| `anchors` | payload 中承载母体连接原子的键；空 = 不收集 |
| `parent_anchor_fields` | `(单, 复)` 锚点字段名，仅供固定 locant 1 的扁平字段 |
| `locant_source` | `attachment` / `attachment_exocyclic` / `anchor_field` |

磺酸族与磺酰胺各一条：

```python
FgSpec("oxoacid", p41=9, path=(0,), anchors=("center_idx",))
FgSpec("sulfonamide", p41=11, path=(1,), anchors=("center_idx",))
```

`layer1/functional_group_inventory.FunctionalGroupClass` 加成员，值必须与 `FgSpec.fg`、`FG_SMARTS` 键逐字一致——`_ANCHOR_KEYS` 与 `_one` 都做 `FunctionalGroupClass(key)`，缺成员即失败。

`FG_ATOM_FNS` 只在特征原子不走「`center_idx` ∪ `surr_idx`」通用形态时登记：`oxoacid` 与 `sulfonamide` 共用 `_oxoacid_atoms`，即「锚点碳 + 中心 P/S + 中心的非碳非氢邻居」。

## 5 P-41 酸式漂移（可选）

该类可整体升为 P-41 类 7 竞争时，用 `fg_registry.oxoacid_is_acid` 判定：中心带 O⁻、或 `oxo_kind in OXO_ACID_KINDS`、或 `n_oh >= OXO_ACID_KIND_BY_H[oxo_kind]`。命中时 `layer2/principal._effective_priority` 换成 `PrincipalPriority(OXO_ACID_P41, path)`，`OXO_ACID_P41 = 8` 排在羧酸之后。三张表同在 `fg_registry`，按新 kind 补行即可。

## 6 母体化接线（L2）

注册完成后 L2 自动接线，只有例外情形才动手：

- `select_principal_group` 取 `PrincipalPriority` 最小的类别作主基团，低优先级 FG 一律成为取代基；互斥由「只选单个最高优先级类」这一结构实现，无需配置互斥键。
- `express_chain_principal` / `express_ring_principal` 按骨架拓扑表达主基团；`_chain_kind` 给出 kind，含氧酸类取 `_oxo_kind_of`，要求该类全部 occurrence 的 `oxo_kind` 一致，不一致即该类不支持。
- `_chain_oxoacid_fields` 写入 `oxo_kind` / `n_oh` / `n_om`；`oxo_kind in _OXO_Z_ANCHORED` 时另做盐门控（金属数与 `n_om` 匹配，门控不过返回 `None` 使候选作废），碳锚定类无门控。
- `_facts` 汇总 `multiplicity` / `relation`（`in_skeleton` / `exocyclic`）/ `attachment_atoms` / `charge_state`。多重度由 occurrence 个数承载，不要另造 kind。
- `assemble` 依次接 hydro 前缀、环内阳离子后缀、`join_kind_name`、阴离子、E/Z、R/S；R/S 的适用集为全部 `FunctionalGroupClass` 值，新 FG 自动纳入。

例外挂钩点：

| 需求 | 挂钩点 |
|---|---|
| 位次由固定 locant 1 的扁平字段承载 | `_SEMANTIC_ANCHOR_FGS` / `_semantic_anchor_fields` 放行并写字段，配 `FgSpec.parent_anchor_fields` |
| L5 整名需要新的计数或盐元数据 | `express_chain_principal` 中按 `kind` 加 producer |
| 该 FG 整体被更高优先级类压制 | `analyzer._SUPPRESSIBLE`（组合 FG）/ `_LEAF_DEMOTED`（降为前缀叶） |

## 7 命名

`layer5/chain_engine._KIND_TABLE` 注册该 kind 的 `_Chain` 条目。字段按源码列全：

```text
词形:   kind / en_suf / zh_suf / coda
位次:   fg / need / no_loc / omit_rule / yl_loc_omit / zh_loc_omit / ene_loc_omit / yne_loc_omit
不饱和: ene_seg / yne_seg / ene_base / yne_suf / unsat_polyol
形态:   ez_ene / ez_ene_multi / wrap / cyclic / cyclic_unsat / aromatic / stem / zh_full
保留名: plain_maps / plain_fn / plain_hook / variant
```

磺酸的条目：

```python
"sulfonic": _Chain(kind="sulfonic", en_suf="sulfonic acid", zh_suf="磺酸", coda="ane",
                   fg="oxoacid", need=1, omit_rule=_omit_term_locant,
                   ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent)
```

- `fg` 是**位次记录的 kind**（= `FgSpec.fg`），不一定是 `kind`：`sulfonic` 的 `fg` 填 `"oxoacid"`，因为 L1 条目键与 L4 记录键都是 `oxoacid`；`sulfonamide` 两者同名。
- `need=1` + `omit_rule=_omit_term_locant` 让单取代短链省位次。
- 多基后缀由 `_generated_mult_fields` 按 `facts.multiplicity` 自动生成，不要为多重度另造 kind。

载荷字段的取值含义：

| 字段 | 取值 / 作用 |
|---|---|
| `coda` | 饱和词干与后缀之间的词尾，`"ane"` 配元音开头的后缀时省 `e` |
| `no_loc` | `"plain"` 出普通名；`"none"` 位次记录缺失即返回 `None` |
| `omit_rule` | `(n, loc, omit) -> bool`，`True` 表示省略位次 |
| `plain_hook` | `(numbered) -> pair`，词尾由分子计数直接决定，不查碳数词表 |
| `variant` | `{scaffold_id: {multiplicity: 覆盖字段}}`，键 `None` 覆盖无 scaffold 情形 |
| `plain_maps` / `plain_fn` | 俗名表 / 派生命名，二者皆空时回落词干拼接 |
| `ene_base` / `yne_suf` | 融合式烯/炔基座；`ene_seg` / `yne_seg` 为段式段形 |
| `unsat_polyol` | 多 FG 词干支持烯/炔插入（多醇/多硫醇） |
| `cyclic` / `cyclic_unsat` | 恒加环前缀 / 仅在不饱和段上加环 |
| `wrap` | `(pair, numbered) -> pair`，整体包裹（E/Z） |

其余词尾落点：

| 场景 | 落点 |
|---|---|
| 中心自任母体的词尾 | `_OXO_TAIL`，键为 `(oxo_kind, 中心酸式氢数)`，由 `plain_hook=_oxoacid_tail` 消费；碳锚定族不查此表 |
| 苯环单取代保留名 | `_BENZENE_RETAINED`（`sulfonic` → benzenesulfonic acid / 苯磺酸） |
| 其他 scaffold 单取代保留名 | `_Chain.variant` 的 `{scaffold_id: {multiplicity: 覆盖字段}}` |
| 环外主基系统名 | `constants.EXO_RING_SUF` 加 `group_class → (单, 复)`；由 `_exo_ring_spec` 改写词尾/词干/位次 |
| 作取代基前缀时的保留名 | `src/namepredict/tools/anchored_table._REGISTRY` 加 `RetainedSubstituent`（如 `sulfo` / `sulfonato`），`anchored` 用带 `*` 的 canonical SMILES |

## 8 原子归属

| 条件 | 落点 |
|---|---|
| 特征原子跨两跳，或中心不是碳骨架成员 | `layer2/parent_select._WHOLE_FG_ATOMS` 加该 FG 类别（现有 `OXOACID`、`SULFONAMIDE`）；命中后 `_kind_fg_atoms` 按 `covered_principal_ids` 整组纳入特征原子 |
| 中心自任母体、臂应退为取代基 | `constants.OXO_CENTER_KINDS` 加该 kind（`phosphate` / `phosphonate` / `sulfate`）；`finalize_parent_ownership` 随即把 `owned_atoms` 取成纯 FG 原子，链仅供编号 |

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
| `_WHOLE_FG_ATOMS` | `layer2/parent_select.py` | 手动（仅跨跳特征原子） |
| `PRINCIPAL_REGISTRY` / kind | `layer2/principal.py`、`principal_expression.py` | 派生 |
| `OXO_CENTER_KINDS` / `EXO_RING_SUF` | `constants.py` | 手动（仅中心母体 / 环外表达） |
| `_KIND_TABLE` / `_OXO_TAIL` / `_BENZENE_RETAINED` | `layer5/chain_engine.py` | 手动 |
| `_FG_LOCANTS` | `layer4/locant_calc.py` | 派生 |
| `_FG_GROUP` | `layer4/locant_calc.py` | 手动（仅位次省略特例） |
| `_REGISTRY` | `src/namepredict/tools/anchored_table.py` | 手动（仅作取代基前缀时） |

## 常见陷阱

1. **漏注册 `FgSpec`**：清单键集与 `_ANCHOR_KEYS` 都来自 `FG_SPECS`，该类在 inventory 中不可见，下游整链静默断掉。
2. **`p41` 填错**：主基团胜负被静默改写，波及已能命名的分子。同类内用 `path` 排先后（`sulfonamide` 的 `path=(1,)` 排在 `amide` 之后）。
3. **特征原子范围决定 L3 提取边界**：范围过宽会把本该作取代基的碳抢走。`_oxoacid_atoms` 只收中心的非碳非氢邻居与锚点碳，碳臂留给链/取代基侧。
4. **`_KIND_TABLE` 与 `_OXO_TAIL` 分离**：前者给该 kind 的基底词形（碳链母体），后者给中心母体按酸式氢数查的词尾。
5. **`_Chain.fg` 不等于 `_Chain.kind`**：`fg` 必须等于 `FgSpec.fg`（位次记录 kind），`sulfonic` 即反例。
6. **环形 / 无环、取代式 / 保留名的分叉要同时覆盖**：保留名走 `variant` 与 `_BENZENE_RETAINED`，环外表达走 `EXO_RING_SUF`，位次省略走 `_FG_GROUP`。
7. **派生表不要手写**：清单键、`_ANCHOR_KEYS`、`PRINCIPAL_REGISTRY`、kind、`_FG_LOCANTS` 都由注册表投影，另立一份会与注册表脱节。

## 验证

1. 简单 FG 作母体：最短链、不同链长各一条，断言词干与位次省略；
2. 不饱和链：验证段式引擎的烯/炔段形态；
3. 与更高优先级 FG 共存：确认本 FG 退为取代基前缀；
4. 与其他取代基叠加：验证编号与位次；
5. 环骨架与环外表达各一条，确认保留名与 `EXO_RING_SUF` 分支被走到。

双语（en / zh）都要断言，用例按主题落 `tests/unit/`。pytest 并行运行；`tests/unit/test_architecture_contracts.py` 校验 `_detect_parts` 键集不越出 `FG_SPECS` 的 `fg` 集。

## 对照示例：中心自任母体（`phosphate`）

与碳锚定路径的三点差异：

- `_oxoacid_entry` 走 `_OXO_Z_ANCHORED` 分支：`center_idx` 即中心原子，臂须单点回接并通过整分子纯度校验，条目不带 `surr_idx`；
- `owned_atoms` 只取 FG 特征原子（`finalize_parent_ownership` 依 `OXO_CENTER_KINDS` 判定），链仅供编号，臂一律退为取代基；
- 词尾走 `plain_hook=_oxoacid_tail` + `_OXO_TAIL`，按 `(oxo_kind, n_oh)` 出 phosphoric acid / dihydrogen phosphate 等，O 侧臂由 `join_kind_name` 按 `ESTER_O_SIDE_KINDS` 拼接。

中心 kind 归入 `constants.OXO_CENTER_KINDS` 后，`namer._apply_salt_suffix` 自动跳过通用盐后缀，盐形态由该类自行组装。
